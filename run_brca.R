library(dndscv)

## With breast cancer dataset

brca_path <- file.path("inst", "extdata", "standard.PCAWG-BRCA-UK_BRCA_27.WGS.tsv")
brca_df <- read.table(brca_path, header = TRUE, sep = "\t")

# Prepare mutation-signatures-data mutation table for dndscv
clean_msd_df <- function(df) {
    filtered_df <- df[df$Mutation.Type == "SBS", ]
    filtered_df <- filtered_df[, c("Sample", "Chromosome", "Start.Position", "Reference.Sequence", "Variant.Sequence")]
    colnames(filtered_df) <- c("sampleID", "chr",      "pos",      "ref",      "alt")
    return(filtered_df)
}

clean_mech_df <- function(df) {
    # Remove extra columns
    filtered_df <- df[, c("sampleID", "chr",      "pos",      "ref",      "alt")]
    return(filtered_df)
}

# Create matrix for 1-parameter substitution model
get_1r_3w_matrix <- function() {
    data(list=sprintf("submod_%s","2r_3w"), package="dndscv")
    mat <- substmodel
    for(i in seq_len(dim(mat)[1])) {
        mat[i, ] <- mat[1, ]
    }
    return(mat)
}

run_subs_models <- function(mut_df) {
    max_coding_muts_per_sample <- 500

    dnds_with_cv_192r = dndscv(mut_df, sm = "192r_3w", max_coding_muts_per_sample = max_coding_muts_per_sample, max_muts_per_gene_per_sample = Inf)
    dnds_with_cv_12r = dndscv(mut_df, sm = "12r_3w", max_coding_muts_per_sample = max_coding_muts_per_sample, max_muts_per_gene_per_sample = Inf)
    dnds_with_cv_2r = dndscv(mut_df, sm = "2r_3w", max_coding_muts_per_sample = max_coding_muts_per_sample, max_muts_per_gene_per_sample = Inf)
    dnds_with_cv_1r = dndscv(mut_df, sm = get_1r_3w_matrix(), max_coding_muts_per_sample = max_coding_muts_per_sample, max_muts_per_gene_per_sample = Inf)
    #dnds_without_cv_192r = dndscv(mut_df, sm = "192r_3w", cv = NULL, max_coding_muts_per_sample = max_coding_muts_per_sample, max_muts_per_gene_per_sample = Inf)
    #dnds_without_cv_12r = dndscv(mut_df, sm = "12r_3w", cv = NULL, max_coding_muts_per_sample = max_coding_muts_per_sample, max_muts_per_gene_per_sample = Inf)
    #dnds_without_cv_2r = dndscv(mut_df, sm = "2r_3w", cv = NULL, max_coding_muts_per_sample = max_coding_muts_per_sample, max_muts_per_gene_per_sample = Inf)
    #dnds_without_cv_1r = dndscv(mut_df, sm = get_1r_3w_matrix(), cv = NULL, max_coding_muts_per_sample = max_coding_muts_per_sample, max_muts_per_gene_per_sample = Inf)

    out <- list(
        "1_with_cv" = dnds_with_cv_1r,
        "2_with_cv" = dnds_with_cv_2r,
        "12_with_cv" = dnds_with_cv_12r,
        "192_with_cv" = dnds_with_cv_192r
        #"1_without_cv" = dnds_without_cv_1r,
        #"2_without_cv" = dnds_without_cv_2r,
        #"12_without_cv" = dnds_without_cv_12r,
        #"192_without_cv" = dnds_without_cv_192r
    )
    return(out)
}

get_sel_df <- function(out, mechanism) {
    out_df <- data.frame()
    for(param_run_name in names(out)) {
        param_run <- out[[param_run_name]]
        param_df <- param_run$sel_cv
        param_df$Num_Rates <- stringr::str_split(param_run_name, "_")[[1]][1]
        param_df$Covariates <- stringr::str_split(param_run_name, "_")[[1]][2]

        out_df <- rbind(out_df, param_df)
    }
    out_df$Mechanism <- mechanism
    out_df$Num_Rates <- factor(out_df$Num_Rates, levels = c("1", "2", "12", "192"))
    return(out_df)
}

get_out_df <- function(out, mechanism) {
    out_df <- data.frame()
    for(param_run_name in names(out)) {
        param_run <- out[[param_run_name]]
        out_df <- rbind(out_df, list(
            Num_Rates = stringr::str_split(param_run_name, "_")[[1]][1],
            AIC = AIC(param_run$poissmodel),
            Covariates = stringr::str_split(param_run_name, "_")[[1]][2],
            globaldnds_all_mle = param_run$globaldnds["wall", "mle"],
            globaldnds_all_cilow = param_run$globaldnds["wall", "cilow"],
            globaldnds_all_cihigh = param_run$globaldnds["wall", "cihigh"]
        ))
    }
    out_df$Num_Rates <- factor(out_df$Num_Rates, levels = c("1", "2", "12", "192"))

    # num_diff_gene_results <- sum(out$`192_with_cv`$sel_cv$qallsubs_cv != out$`192_without_cv`$sel_cv$qallsubs_cv)
    # if(num_diff_gene_results > 0) {
    #     warning("Results are different with and without covariates!!")
    # } else {
    #     # FOR NOW, only use without covariates since the results seem to be the same
    #     out_df <- out_df[out_df$Covariates == "without"]
    # }
    out_df$Mechanism <- mechanism
    return(out_df)
}

brca_clean_df <- clean_msd_df(brca_df)

out <- run_subs_models(brca_clean_df)
saveRDS(out, "data/out.rds")


## Change transition mutations to transversions

add_mechanism_cols <- function(df) {
    df$mechanism <- apply(df, 1, function(row) {
        mut_ref <- row[["ref"]]
        mut_var <- row[["alt"]]
        transition <- c( "AG", "GA", "CT", "TC" )
        transversion <- c( "AT", "TA", "AC", "CA", "CG", "GC", "TG", "GT" )
        if(paste0(mut_ref, mut_var) %in% transition) {
            return("Transition")
        } else if(paste0(mut_ref, mut_var) %in% transversion) {
            return("Transversion")
        } else {
            stop("Error, neither transition nor transversion")
        }
    })
    df$is_transition <- apply(df, 1, function(row) {
        return(row[["mechanism"]] == "Transition")
    })
    return(df)
}

# Ref mapping to any var that would make it a transversion
ts_to_tv_map <- list(
    `A` = c("T", "C"),
    `C` = c("A", "G"),
    `G` = c("C", "T"),
    `T` = c("G", "A")
)
# Ref mapping to any var that would make it a transition
tv_to_ts_map <- list(
    `A` = "G",
    `C` = "T",
    `G` = "A",
    `T` = "C"
)

convert_ts_to_tv <- function(df) {
    df <- data.frame(df)
    df$alt <- apply(df, 1, function(row) {
        mut_ref <- row[["ref"]]
        mut_var <- row[["alt"]]
        is_transition <- as.logical(row[["is_transition"]])
        if(!is_transition) {
            return(mut_var)
        } else {
            # This is a transition, convert to transversion.
            # Randomly select one of the two options for the given ref base.
            return(sample(ts_to_tv_map[[mut_ref]], 1))
        }
    })
    df <- add_mechanism_cols(df)
    return(df)
}
convert_tv_to_ts <- function(df) {
    df <- data.frame(df)
    df$alt <- apply(df, 1, function(row) {
        mut_ref <- row[["ref"]]
        mut_var <- row[["alt"]]
        is_transition <- as.logical(row[["is_transition"]])
        if(!is_transition) {
            # This is a transversion, convert to transition.
            # Only one choice per ref base.
            return(tv_to_ts_map[[mut_ref]])
        } else {
            return(mut_var)
        }
    })
    df <- add_mechanism_cols(df)
    return(df)
}

brca_clean_df <- add_mechanism_cols(brca_clean_df)

brca_transversion_only_df <- convert_ts_to_tv(brca_clean_df)
brca_transition_only_df <- convert_tv_to_ts(brca_clean_df)

brca_transversion_only_clean_df <- clean_mech_df(brca_transversion_only_df)
brca_transition_only_clean_df <- clean_mech_df(brca_transition_only_df)

out_transversion_only <- run_subs_models(brca_transversion_only_clean_df)
saveRDS(out_transversion_only, "data/out_transversion_only.rds")

out_transition_only <- run_subs_models(brca_transition_only_clean_df)
saveRDS(out_transition_only, "data/out_transition_only.rds")


