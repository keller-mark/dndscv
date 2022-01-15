
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

transition_missense_only_df <- transition_only_df[transition_only_df$impact %in% c("Synonymous", "Missense"), ]
transition_nonsense_only_df <- transition_only_df[transition_only_df$impact %in% c("Synonymous", "Nonsense"), ]
transversion_missense_only_df <- transversion_only_df[transversion_only_df$impact %in% c("Synonymous", "Missense"), ]
transversion_nonsense_only_df <- transversion_only_df[transversion_only_df$impact %in% c("Synonymous", "Nonsense"), ]

# Count dN/dS ratio
get_dN <- function(df, mechanism_types, ns_types) {
    ns_df <- df[df$impact %in% ns_types, ]
    ns_subs <- sum(!ns_df[ns_df$mechanism %in% mechanism_types, ]$is_synonymous)
    ns_sites <- sum(!ns_df$is_synonymous)
    return(ns_subs / ns_sites)
}
get_dS <- function(df, mechanism_types, ns_types) {
    s_subs <- sum(df[df$mechanism %in% mechanism_types, ]$is_synonymous)
    s_sites <- sum(df$is_synonymous)
    return(s_subs / s_sites)
}
get_dNdS <- function(df, mechanism_types, ns_types) {
    return(get_dN(df, mechanism_types, ns_types) / get_dS(df, mechanism_types, ns_types))
}

# Plot
library(ggplot2)

dNdS_df <- data.frame(
    Mechanism = c(
        "All",
        "Transition Only",
        "Transition Only (Missense)",
        "Transition Only (Nonsense)",
        "Transversion Only",
        "Transversion Only (Missense)",
        "Transversion Only (Nonsense)"
    ),
    dNdS = c(
        get_dNdS(mut_base_df, c("Transition", "Transversion"), c("Missense", "Nonsense")),
        get_dNdS(mut_base_df, c("Transition"), c("Missense", "Nonsense")),
        get_dNdS(mut_base_df, c("Transition"), c("Missense")),
        get_dNdS(mut_base_df, c("Transition"), c("Nonsense")),
        get_dNdS(mut_base_df, c("Transversion"), c("Missense", "Nonsense")),
        get_dNdS(mut_base_df, c("Transversion"), c("Missense")),
        get_dNdS(mut_base_df, c("Transversion"), c("Nonsense"))
    ),
    Num_NS = c(
        get_dN(mut_base_df, c("Transition", "Transversion"), c("Missense", "Nonsense")),
        get_dN(mut_base_df, c("Transition"), c("Missense", "Nonsense")),
        get_dN(mut_base_df, c("Transition"), c("Missense")),
        get_dN(mut_base_df, c("Transition"), c("Nonsense")),
        get_dN(mut_base_df, c("Transversion"), c("Missense", "Nonsense")),
        get_dN(mut_base_df, c("Transversion"), c("Missense")),
        get_dN(mut_base_df, c("Transversion"), c("Nonsense"))
    ),
    Num_S = c(
        get_dS(mut_base_df, c("Transition", "Transversion"), c("Missense", "Nonsense")),
        get_dS(mut_base_df, c("Transition"), c("Missense", "Nonsense")),
        get_dS(mut_base_df, c("Transition"), c("Missense")),
        get_dS(mut_base_df, c("Transition"), c("Nonsense")),
        get_dS(mut_base_df, c("Transversion"), c("Missense", "Nonsense")),
        get_dS(mut_base_df, c("Transversion"), c("Missense")),
        get_dS(mut_base_df, c("Transversion"), c("Nonsense"))
    )
)
dNdS_df$Ratio = paste(round(dNdS_df$Num_NS, 2), round(dNdS_df$Num_S, 2), sep = " / ")
dNdS_df$Mechanism <- factor(dNdS_df$Mechanism, levels = dNdS_df$Mechanism)

ggplot(dNdS_df, aes(Mechanism, dNdS, color = Mechanism, label = Ratio)) +
    geom_point(size = 4) +
    geom_text(check_overlap = TRUE, nudge_y = 0.1, nudge_x = 0) +
    ylim(0, 2) +
    geom_hline(yintercept=1.0, linetype="dashed", color = "black") +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1))

ggsave("pqe_plots/dNdS_no_control_simplest.png", scale = 1, width = 6, height = 6, units = "in")

# TODO: incorporate a substitution weight matrix


