# Tutorial
# Reference: http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html

library(dndscv)
library(ggplot2)

# Comparing results using different substitution models

data("dataset_normalskin", package="dndscv")
data("dataset_normalskin_genes", package="dndscv")

compare_subs_models <- function(mut_df) {
    return(comparison_df)
}
mut_df <- m
# Skin dataset analysis with full 192-category model
dnds_192r = dndscv(mut_df, sm = "192r_3w")
# Analysis with 12-category model
dnds_12r = dndscv(mut_df, sm = "12r_3w")
# Skin dataset analysis with transition-transversion model
dnds_2r = dndscv(mut_df, sm = "2r_3w")

comparison_df <- data.frame(
    Substitution_Model = c(
        "HKY85",
        "GTR",
        "192"
    ),
    AIC = c(
        AIC(dnds_2r$poissmodel),
        AIC(dnds_12r$poissmodel),
        AIC(dnds_192r$poissmodel)
    ),
    wall_mle = c(
        dnds_2r$globaldnds["wall", "mle"],
        dnds_12r$globaldnds["wall", "mle"],
        dnds_192r$globaldnds["wall", "mle"]
    ),
    wall_cilow = c(
        dnds_2r$globaldnds["wall", "cilow"],
        dnds_12r$globaldnds["wall", "cilow"],
        dnds_192r$globaldnds["wall", "cilow"]
    ),
    wall_cihigh = c(
        dnds_2r$globaldnds["wall", "cihigh"],
        dnds_12r$globaldnds["wall", "cihigh"],
        dnds_192r$globaldnds["wall", "cihigh"]
    )
)
comparison_df$Substitution_Model <- factor(comparison_df$Substitution_Model, levels = c("HKY85", "GTR", "192"))

ggplot(comparison_df, aes(Substitution_Model, AIC, color = Substitution_Model)) +
    geom_point(size = 4) +
    theme(legend.position="none") +
    ggtitle("AIC by substitution model on skin dataset")

ggsave("pqe_plots/AIC_by_sm_skin.png", scale = 1, width = 4, height = 4, units = "in")
