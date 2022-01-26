library(parallel)

sigs_path <- "inst/extdata/COSMIC-signatures.SBS-96.tsv"
sigs_df <- t(read.table(sigs_path, header = TRUE, row.names = 1, check.names = FALSE))

sbs1 <- as.list(sigs_df[, "1"])
sbs22 <- as.list(sigs_df[, "22"])

N_SAMPLE <- 50
N_MUT_PER_SAMPLE <- 50

# SETUP

# https://github.com/pkerpedjiev/negspy/blob/master/negspy/data/hg19/chromInfo.txt
chr_len <- c(
    chr1 =	249250621,
    chr2 =	243199373,
    chr3 =	198022430,
    chr4 =	191154276,
    chr5 =	180915260,
    chr6 =	171115067,
    chr7 =	159138663,
    chr8 =	146364022,
    chr9 =	141213431,
    chr10 =	135534747,
    chr11 =	135006516,
    chr12 =	133851895,
    chr13 =	115169878,
    chr14 =	107349540,
    chr15 =	102531392,
    chr16 =	90354753,
    chr17 =	81195210,
    chr18 =	78077248,
    chr19 =	59128983,
    chr20 =	63025520,
    chr21 =	48129895,
    chr22 =	51304566,
    chrX =	155270560,
    chrY =	59373566
)
genome_len <- sum(chr_len)

data("refcds_hg19", package="dndscv")
gene_list <- sapply(RefCDS, function(x) x$gene_name)

# Expanding the reference sequences [for faster access]
for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$seq_cds <- base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
    RefCDS[[j]]$seq_cds1up <- base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
    RefCDS[[j]]$seq_cds1down <- base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
    if (!is.null(RefCDS[[j]]$seq_splice)) {
        RefCDS[[j]]$seq_splice <- base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
        RefCDS[[j]]$seq_splice1up <- base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
        RefCDS[[j]]$seq_splice1down <- base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
    }
}
ind <- setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
gr_genes_ind <- ind[gr_genes$names]

# Subfunction: obtaining the codon positions of a coding mutation given the exon intervals
chr2cds <- function(pos,cds_int,strand) {
    if (strand==1) {
        return(which(unlist(apply(cds_int, 1, function(x) x[1]:x[2])) %in% pos))
    } else if (strand==-1) {
        return(which(rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2]))) %in% pos))
    }
}

## END SETUP


find_next_match <- function(ref_trinuc, mut_abs_pos) {
    # Does this position fall in a coding region?

    prev_chroms <- !(mut_abs_pos <= cumsum(chr_len))
    mut_chr_ind <- sum(prev_chroms) + 1
    mut_chr <- names(chr_len)[mut_chr_ind]
    mut_pos <- mut_abs_pos - sum(chr_len[prev_chroms])

    mut_chr_hg19 <- substr(mut_chr, 4, stringr::str_length(mut_chr))

    # Mapping mutations to genes
    gr_muts <- GenomicRanges::GRanges(c(mut_chr_hg19), IRanges::IRanges(c(mut_pos),c(mut_pos)))
    ol <- as.data.frame(GenomicRanges::findOverlaps(gr_muts, gr_genes, type="any", select="all"))

    # If no overlap with a coding region, stop.
    if(dim(ol)[1] == 0) {
        return(FALSE)
    }

    ol_gene_ind <- gr_genes_ind[ol[,2]]
    ol_gene_info <- RefCDS[ol_gene_ind][[1]]

    if (any(mut_pos == ol_gene_info$intervals_splice)) { # Essential splice-site substitution
        return(FALSE)
    } else { # Coding substitution
        count <- 0
        repeat {
            pos_ind <- chr2cds(mut_pos, ol_gene_info$intervals_cds, ol_gene_info$strand)
            if(length(pos_ind) == 0) {
                return(FALSE)
            }
            ref3 <- sprintf("%s%s%s", ol_gene_info$seq_cds1up[pos_ind], ol_gene_info$seq_cds[pos_ind], ol_gene_info$seq_cds1down[pos_ind])

            if(ref3 == ref_trinuc) {
                return(list(chr = mut_chr_hg19, pos = mut_pos, strand = ol_gene_info$strand))
            }

            mut_pos <- mut_pos + 1

            if(count > 10000) {
                return(FALSE)
            }
            count <- count + 1
        }
    }
}

make_sample_mut_df <- function(sample_id, n_mut, mut_sig) {
    df <- data.frame()
    for(mut_i in seq_len(n_mut)) {
        ref_cat <- sample(names(mut_sig), 1, prob = as.numeric(mut_sig))
        ref_cat_vec <- strsplit(ref_cat, "")[[1]]
        ref_trinuc <- sprintf("%s%s%s", ref_cat_vec[1], ref_cat_vec[3], ref_cat_vec[7])
        ref <- ref_cat_vec[3]
        alt <- ref_cat_vec[5]
        next_match <- NULL
        repeat {
            mut_abs_pos <- floor(runif(n = 1, min = 1, max = genome_len))
            next_match <- find_next_match(ref_trinuc, mut_abs_pos)
            if(is.list(next_match)) {
                break
            }
        }
        df <- rbind(df, list(
            pos = next_match$pos,
            chr = next_match$chr,
            strand = next_match$strand,
            alt = alt,
            ref = ref
        ))
        print(mut_i)
    }

    df$sampleID <- sample_id
    return(df)
}

make_sbs1_mut_df_for_sample <- function(sample_id) {
    return(make_sample_mut_df(sample_id, N_MUT_PER_SAMPLE, sbs1))
}
make_sbs22_mut_df_for_sample <- function(sample_id) {
    return(make_sample_mut_df(sample_id, N_MUT_PER_SAMPLE, sbs22))
}

# Make output dataframes

sbs1_samples <- as.list(paste0("sbs1_sample_v2_", seq_len(N_SAMPLE)))
sbs22_samples <- as.list(paste0("sbs22_sample_v2_", seq_len(N_SAMPLE)))

sbs1_mut_dfs <- mclapply(sbs1_samples, make_sbs1_mut_df_for_sample, mc.cores = detectCores())
sbs1_mut_df <- do.call(rbind, sbs1_mut_dfs)
write.csv(sbs1_mut_df, file = "inst/extdata/neutral_sbs1_v2.csv")

sbs22_mut_dfs <- mclapply(sbs22_samples, make_sbs22_mut_df_for_sample, mc.cores = detectCores())
sbs22_mut_df <- do.call(rbind, sbs22_mut_dfs)
write.csv(sbs22_mut_df, file = "inst/extdata/neutral_sbs22_v2.csv")

