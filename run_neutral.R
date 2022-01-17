library(dndscv)

## Run on neutral datasets

dna_complement <- list(
    `A` = "T",
    `C` = "G",
    `G` = "C",
    `T` = "A"
)

# Prepare table for dndscv
clean_sbs_df <- function(df) {
    # Take reverse complement if on -1 strand
    df$ref <- apply(df, 1, function(row) {
        if(as.numeric(row[["strand"]]) == -1) {
            return(dna_complement[[row[["ref"]]]])
        } else {
            return(row[["ref"]])
        }
    })
    df$alt <- apply(df, 1, function(row) {
        if(as.numeric(row[["strand"]]) == -1) {
            return(dna_complement[[row[["alt"]]]])
        } else {
            return(row[["alt"]])
        }
    })

    # Remove extra columns
    df <- df[, c("sampleID", "chr",      "pos",      "ref",      "alt")]

    return(df)
}

sbs1_path <- "inst/extdata/neutral_sbs1_v2.csv"
sbs1_df <- clean_sbs_df(read.table(sbs1_path, header = TRUE, sep = ","))

sbs22_path <- "inst/extdata/neutral_sbs22_v2.csv"
sbs22_df <- clean_sbs_df(read.table(sbs22_path, header = TRUE, sep = ","))

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

out_sbs1 <- run_subs_models(sbs1_df)
saveRDS(out_sbs1, "inst/extdata/out_sbs1.rds")

rm("out_sbs1")

out_sbs22 <- run_subs_models(sbs22_df)
saveRDS(out_sbs22, "inst/extdata/out_sbs22.rds")

