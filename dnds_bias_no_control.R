# Simulation without controlling for trinucleotide context

dom_sig <- list(
    "TCG" = c(
        "A" = 1.0,
        "G" = 0.0,
        "T" = 0.0
    ),

)

sigs_path <- file.path("inst", "extdata", "COSMIC-signatures.SBS-96.tsv")
sigs_df <- t(read.table(sigs_path))
dom_sig <- sigs_df[]

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
gene_list = sapply(RefCDS, function(x) x$gene_name)

# Expanding the reference sequences [for faster access]
for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
    RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
    RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
    if (!is.null(RefCDS[[j]]$seq_splice)) {
        RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
        RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
        RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
    }
}
ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
gr_genes_ind = ind[gr_genes$names]

# 1. Choose genomic location
# Absolute genomic position for mutation
mut_abs_pos <- runif(n = 1, min = 1, max = genome_len)
# TODO : delete next lines
mut_abs_pos <- 249250621 + 243199373 + 2 # should be on chr3

mut_abs_pos <- sum(chr_len[1:16]) + 7572927 + 10 # TP53, should be on chr17

prev_chroms <- !(mut_abs_pos <= cumsum(chr_len))
mut_chr_ind <- sum(prev_chroms) + 1
mut_chr <- names(chr_len)[mut_chr_ind]
mut_pos <- mut_abs_pos - sum(chr_len[prev_chroms])

# Mapping mutations to genes
gr_muts = GenomicRanges::GRanges(c(mut_chr_ind), IRanges::IRanges(c(mut_pos),c(mut_pos)))
ol = as.data.frame(GenomicRanges::findOverlaps(gr_muts, gr_genes, type="any", select="all"))

# TODO: if no overlap with a coding region, skip
ol_gene_ind <- gr_genes_ind[ol[,2]]
ol_gene_info <- RefCDS[ol_gene_ind]
ol_gene_range <- gr_genes[ol[,2]]

# TODO: consider strand?, could have an overlap from each strand or more than one
ol_start <- GenomicRanges::start(ol_gene_range[1])
mut_gene_pos <- mut_pos - ol_start # position of the mutation from the start of the gene

# Compute impact of mutation
geneind = ol_gene_ind
pos = mut_pos

# Subfunction: obtaining the codon positions of a coding mutation given the exon intervals
chr2cds = function(pos,cds_int,strand) {
    if (strand==1) {
        return(which(unlist(apply(cds_int, 1, function(x) x[1]:x[2])) %in% pos))
    } else if (strand==-1) {
        return(which(rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2]))) %in% pos))
    }
}

# Annotating the functional impact of each substitution and populating the N matrices
impind = NA
ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = codonsub = array(NA, 1)
j = 1

if (any(pos == RefCDS[[geneind]]$intervals_splice)) { # Essential splice-site substitution

    impact[j] = "Essential_Splice"; impind = 4
    pos_ind = (pos==RefCDS[[geneind]]$intervals_splice)
    cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
    ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind])
    mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
    aachange[j] = ntchange[j] = codonsub[j] = "."

} else { # Coding substitution

    pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
    cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
    ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
    mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
    codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
    old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
    pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
    new_codon = old_codon; new_codon[pos_in_codon] = mutations$mut_cod[j]
    old_aa = seqinr::translate(old_codon, numcode = numcode)
    new_aa = seqinr::translate(new_codon, numcode = numcode)
    aachange[j] = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
    ntchange[j] = sprintf('%s%0.0f%s',mutations$ref_cod[j],pos_ind,mutations$mut_cod[j])
    codonsub[j] = sprintf('%s>%s',paste(old_codon,collapse=""),paste(new_codon,collapse=""))

    # Annotating the impact of the mutation
    if (new_aa == old_aa){
        impact[j] = "Synonymous"; impind = 1
    } else if (new_aa == "*"){
        impact[j] = "Nonsense"; impind = 3
    } else if (old_aa != "*"){
        impact[j] = "Missense"; impind = 2
    } else if (old_aa=="*") {
        impact[j] = "Stop_loss"; impind = NA
    }
}


