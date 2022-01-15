# Tutorial
# Reference: http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html

library(dndscv)

data("dataset_simbreast", package="dndscv")
dndsout = dndscv(mutations)
print(dndsout$globaldnds)

# Using different substitution models

# 192 rates (used as default)
data("submod_192r_3w", package="dndscv")
colnames(substmodel) = c("syn","mis","non","spl")
head(substmodel)
dim(substmodel)

# 12 rates (no context-dependence)
data("submod_12r_3w", package="dndscv")
colnames(substmodel) = c("syn","mis","non","spl")
head(substmodel)

# 2 rates (classic ts/tv model)
data("submod_2r_3w", package="dndscv")
colnames(substmodel) = c("syn","mis","non","spl")
head(substmodel)

# Using dndscv with a different substitution model

data("dataset_normalskin", package="dndscv")
data("dataset_normalskin_genes", package="dndscv")
dndsskin_2r = dndscv(m, gene_list=target_genes, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, sm = "2r_3w")

print(dndsskin_2r$mle_submodel)





