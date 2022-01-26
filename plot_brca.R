# Plot results

library(ggplot2)
library(ggrepel)

## With breast cancer dataset
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

# Load from cached
out <- readRDS("inst/extdata/out.rds")
out_df <- get_out_df(out, "Original")

ggplot(out_df, aes(Num_Rates, AIC, color = Num_Rates)) +
    geom_point(size = 4) +
    theme(legend.position="none") +
    ggtitle("AIC by substitution model on PCAWG-BRCA-UK dataset")

ggsave("pqe_plots/AIC_by_sm_brca_uk.png", scale = 1, width = 6, height = 5, units = "in")

ggplot(out_df, aes(x = Num_Rates, y = globaldnds_all_mle, group = Num_Rates)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin=globaldnds_all_cilow, ymax=globaldnds_all_cihigh), width=.2)  +
    expand_limits(y = 0) +
    geom_hline(yintercept = 1.0, linetype="dashed", color = "black") +
    theme(legend.position="none") +
    labs(
        title="Global dN/dS by substitution model on PCAWG-BRCA-UK dataset",
        x="Substitution Model: Number of Rate Parameters",
        y = "Global dN/dS",
        subtitle = "Original dataset"
    )

ggsave("pqe_plots/dnds_by_sm_brca_uk.png", scale = 1, width = 6, height = 5, units = "in")

## Change transition mutations to transversions

out_transversion_only <- readRDS("inst/extdata/out_transversion_only.rds")
out_transition_only <- readRDS("inst/extdata/out_transition_only.rds")

# Process cached
out_transversion_only_df <- get_out_df(out_transversion_only, "Transversion Only")
out_transition_only_df <- get_out_df(out_transition_only, "Transition Only")

out_all_df <- rbind(out_transversion_only_df, out_transition_only_df)


## Plot

ggplot(out_all_df, aes(x = Num_Rates, y = globaldnds_all_mle, color = Mechanism)) +
    geom_point(size = 3, aes(color = Mechanism), position = position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin=globaldnds_all_cilow, ymax=globaldnds_all_cihigh, color = Mechanism), width=.2, position = position_dodge(width=0.5))  +
    expand_limits(y = 0) +
    geom_hline(yintercept = 1.0, linetype="dashed", color = "black") +
    labs(
        title="Global dN/dS by substitution model on PCAWG-BRCA-UK dataset",
        x="Substitution Model: Number of Rate Parameters",
        y = "Global dN/dS",
        color = "Dataset"
    )

ggsave("pqe_plots/dnds_by_sm_tstv_brca_uk.png", scale = 1, width = 8, height = 6, units = "in")


# Plot gene-level data
sel_df <- get_sel_df(out, "Original")
sel_transversion_only_df <- get_sel_df(out_transversion_only, "Transversion Only")
sel_transition_only_df <- get_sel_df(out_transition_only, "Transition Only")


sel_all_df <- rbind(sel_transversion_only_df, sel_transition_only_df)

ggplot(sel_all_df, aes(x = Num_Rates, y = qallsubs_cv, color = Mechanism)) +
    geom_boxplot() +
    ylim(0, 1) +
    labs(title="Distribution of gene-level significance values", y = "q-value", x = "Substitution Model: Number of Rate Parameters", color = "Dataset")

# Volcano plots
sel_all_df$minuslog10p_all <- apply(sel_all_df, 1, function(row) {
    return(-log10(as.numeric(row[["qallsubs_cv"]])))
})
sel_all_df$minuslog10p_mis <- apply(sel_all_df, 1, function(row) {
    return(-log10(as.numeric(row[["qmis_cv"]])))
})
sel_all_df$minuslog10p_trunc <- apply(sel_all_df, 1, function(row) {
    return(-log10(as.numeric(row[["qtrunc_cv"]])))
})

mis_signif_df <- sel_all_df[sel_all_df$qallsubs_cv < 0.1 & sel_all_df$minuslog10p_mis != Inf & sel_all_df$wmis_cv != 0.0, ]
trunc_signif_df <- sel_all_df[sel_all_df$qallsubs_cv < 0.1 & sel_all_df$minuslog10p_trunc != Inf & sel_all_df$wnon_cv != 0.0, ]

ggplot(mis_signif_df, aes(x = log2(wmis_cv), y = minuslog10p_mis, color = Mechanism, label = gene_name)) +
    geom_point(aes(shape = Num_Rates), size = 2) +
    geom_text_repel(show.legend = FALSE) +
    labs(
        title="Distribution of gene-level significance values",
        y = "-log10(missense q-value)",
        x = "log2(missense dN/dS)",
        color = "Dataset"
    )

ggsave("pqe_plots/dnds_by_sm_tstv_brca_uk_volcano_mis.png", scale = 1, width = 8, height = 6, units = "in")

ggplot(trunc_signif_df, aes(x = log2(wnon_cv), y = minuslog10p_trunc, color = Mechanism, label = gene_name)) +
    geom_point(aes(shape = Num_Rates), size = 2) +
    geom_text_repel(show.legend = FALSE) +
    labs(
        title="Distribution of gene-level significance values",
        y = "-log10(truncating q-value)",
        x = "log2(nonsense dN/dS)",
        color = "Dataset"
    )

ggsave("pqe_plots/dnds_by_sm_tstv_brca_uk_volcano_non.png", scale = 1, width = 8, height = 6, units = "in")
