# Plot results

library(ggplot2)
library(ggrepel)


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
out_sbs1 <- readRDS("inst/extdata/out_sbs1.rds")
out_sbs22 <- readRDS("inst/extdata/out_sbs22.rds")

out_sbs1_df <- get_out_df(out_sbs1, "SBS1")
out_sbs22_df <- get_out_df(out_sbs22, "SBS22")

out_all_df <- rbind(out_sbs1_df, out_sbs22_df)

## Plot

ggplot(out_all_df, aes(x = Num_Rates, y = globaldnds_all_mle, color = Mechanism)) +
    geom_point(size = 3, aes(color = Mechanism), position = position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin=globaldnds_all_cilow, ymax=globaldnds_all_cihigh, color = Mechanism), width=.2, position = position_dodge(width=0.5))  +
    expand_limits(y = 0) +
    geom_hline(yintercept = 1.0, linetype="dashed", color = "black") +
    labs(
        title="Global dN/dS by substitution model on simulated dataset",
        x="Substitution Model: Number of Rate Parameters",
        y = "Global dN/dS",
        color = "Dataset"
    )

ggsave("pqe_plots/dnds_by_sm_tstv_using_sigs.png", scale = 1, width = 6, height = 5, units = "in")


# Plot gene-level data
sel_sbs1_df <- get_sel_df(out_sbs1, "SBS1")
sel_sbs22_df <- get_sel_df(out_sbs22, "SBS22")

sel_all_df <- rbind(sel_sbs1_df, sel_sbs22_df)

# None is significant
