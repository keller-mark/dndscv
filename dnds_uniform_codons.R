
codon_path <- "inst/extdata/codon.txt"
codon_df <- read.table(codon_path)
colnames(codon_df) <- c("Codon",	    "AminoAcid",	"Letter",	 "FullName")

codon_df$Base1 <- lapply(codon_df$Codon, function(x) { substr(x, 1, 1) })
codon_df$Base2 <- lapply(codon_df$Codon, function(x) { substr(x, 2, 2) })
codon_df$Base3 <- lapply(codon_df$Codon, function(x) { substr(x, 3, 3) })

cross_df_with_vec <- function(orig_df, new_colname, cross_vec) {
    out_df <- NA

    for(val in cross_vec) {
        df_copy <- data.frame(orig_df)
        df_copy[[new_colname]] <- val
        if(is.na(out_df)) {
            out_df <- df_copy
        } else {
            out_df <- rbind(out_df, df_copy)
        }
    }

    out_df
}

mut_pos_df <- cross_df_with_vec(codon_df, "mut_pos", c(1, 2, 3))
mut_base_df <- cross_df_with_vec(mut_pos_df, "mut_base_idx", c(1, 2, 3))

base_cols <- c("Base1", "Base2", "Base3")
nt <- c("A", "C", "G", "T")

mut_base_df$mut_ref <- apply(mut_base_df, 1, function(row) {
    return(row[[base_cols[row$mut_pos]]])
})
mut_base_df$mut_var <- apply(mut_base_df, 1, function(row) {
    remaining_nt <- setdiff(nt, c(row$mut_ref))
    return(remaining_nt[row$mut_base_idx])
})
mut_base_df$mut_codon <- apply(mut_base_df, 1, function(row) {
    ref_codon <- strsplit(row$Codon, "")[[1]]
    mut_pos <- row$mut_pos
    ref_codon[mut_pos] <- row$mut_var
    return(paste(ref_codon, collapse = ""))
})
mut_base_df$impact <- apply(mut_base_df, 1, function(row) {
    old_codon <- row$Codon
    new_codon <- row$mut_codon
    old_aa = seqinr::translate(seqinr::s2c(old_codon), numcode = 1)
    new_aa = seqinr::translate(seqinr::s2c(new_codon), numcode = 1)
    # Annotating the impact of the mutation
    if (new_aa == old_aa){
        impact <- "Synonymous";
    } else if (new_aa == "*"){
        impact <- "Nonsense";
    } else if (old_aa != "*"){
        impact <- "Missense";
    } else if (old_aa=="*") {
        impact <- "Stop_loss";
    }
})
mut_base_df$is_synonymous <- apply(mut_base_df, 1, function(row) {
    return(row$impact == "Synonymous")
})
mut_base_df$mechanism <- apply(mut_base_df, 1, function(row) {
    mut_ref <- row$mut_ref
    mut_var <- row$mut_var
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
mut_base_df$is_transition <- apply(mut_base_df, 1, function(row) {
    return(row$mechanism == "Transition")
})

# Count number of ways to make a transition
transition_only_df <- mut_base_df[mut_base_df$is_transition, ]
transversion_only_df <- mut_base_df[!mut_base_df$is_transition, ]

# Count dN/dS ratio
dNdS_overall <- sum(!mut_base_df$is_synonymous) / sum(mut_base_df$is_synonymous)
dNdS_transition_only <- sum(!transition_only_df$is_synonymous) / sum(transition_only_df$is_synonymous)
dNdS_transversion_only <- sum(!transversion_only_df$is_synonymous) / sum(transversion_only_df$is_synonymous)

