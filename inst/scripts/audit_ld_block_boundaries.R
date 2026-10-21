#!/usr/bin/env Rscript
# Worker for the LD block-boundary audit: computes genotype LD for one
# superpopulation in one sliding window and saves every HapMap3 SNP pair within
# `step_bp` of each other, so summarize_ld_block_audit.R can measure how much LD
# any candidate block partition discards at its boundaries.
#
# Why a separate genotype computation: the HDF5 panels only store within-block
# LD, so they cannot show what a boundary throws away. LD here comes from
# computeSuperPopLDWindow() -- the same corrected estimator the panels use --
# and is then restricted to PRS-CS's HapMap3 SNP set, the SNPs PRS-CS actually
# analyzes.
#
# Windows are [start0 + k*step, start0 + (k+2)*step - 1], with start0 = the
# chromosome's first HapMap3 bp rounded down to a multiple of step. Every pair
# (a <= b) with b - a <= step lies inside the window whose first half contains
# a, so each pair is kept only in that window -- no pair is counted twice.
#
# Usage:
#   Rscript audit_ld_block_boundaries.R <33kg_index.gz> <33kg_geno.gz> \
#     <33kg_pop_desc.txt> <hm3_snpinfo> <chr> <superpop> <window_k> <out_rds> \
#     [step_bp=1000000] [maf_cutoff=0.01]
#
#   <hm3_snpinfo>  PRS-CS's snpinfo_1kg_hm3 (CHR SNP BP A1 A2 MAF, hg19).
#   <window_k>     0-based window index (the SLURM array task id).
#
# Output: an RDS list with `pairs` (data frame i, j, r; i/j index rows of the
# chromosome's HapMap3 table in bp order), `present` (HapMap3 row indices that
# passed QC and lie in this window's first half), `n_samples`, and metadata.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript audit_ld_block_boundaries.R <index> <geno> <pop_desc> <hm3_snpinfo> <chr> <superpop> <window_k> <out_rds> [step_bp=1000000] [maf_cutoff=0.01]")
}
index_file <- args[[1]]
geno_file <- args[[2]]
pop_desc_file <- args[[3]]
hm3_file <- args[[4]]
chr <- as.integer(args[[5]])
superpop <- toupper(args[[6]])
k <- as.integer(args[[7]])
out_rds <- args[[8]]
step_bp <- if (length(args) >= 9) as.numeric(args[[9]]) else 1e6
maf_cutoff <- if (length(args) >= 10) as.numeric(args[[10]]) else 0.01

suppressPackageStartupMessages(library(gauss))

hm3 <- utils::read.table(hm3_file, header = TRUE, stringsAsFactors = FALSE)
hm3 <- hm3[hm3$CHR == chr, , drop = FALSE]
hm3 <- hm3[order(hm3$BP), , drop = FALSE]
rownames(hm3) <- NULL

start0 <- floor(min(hm3$BP) / step_bp) * step_bp
win_start <- start0 + k * step_bp
win_end <- win_start + 2 * step_bp - 1
anchor_end <- win_start + step_bp - 1

pop_desc <- utils::read.table(pop_desc_file, header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, quote = "")
n_samples <- sum(pop_desc[[2]][toupper(pop_desc[[3]]) == superpop])

empty <- list(
  pairs = data.frame(i = integer(0), j = integer(0), r = numeric(0)),
  present = integer(0), n_samples = n_samples, chr = chr, superpop = superpop,
  window_k = k, win_start = win_start, win_end = win_end, step_bp = step_bp
)

message(sprintf("chr%d %s window %d: [%s, %s]", chr, superpop, k,
  format(win_start, big.mark = ","), format(win_end, big.mark = ",")))

t0 <- Sys.time()
res <- tryCatch(
  computeSuperPopLDWindow(
    chr = chr, start_bp = win_start, end_bp = win_end,
    superpopulation_label = superpop,
    reference_index_file = index_file,
    reference_data_file = geno_file,
    reference_pop_desc_file = pop_desc_file,
    maf_cutoff = maf_cutoff
  ),
  error = function(e) {
    # Gaps (e.g. the centromere) leave too few SNPs to compute LD.
    if (grepl("Fewer than 2 SNPs", conditionMessage(e))) NULL else stop(e)
  }
)
if (is.null(res)) {
  message("Too few SNPs in window; writing an empty result.")
  saveRDS(empty, out_rds)
  quit(save = "no")
}

# Match to HapMap3 by position and allele pair (either orientation); rsids in
# the reference index can be placeholders like ".".
snp <- res$snplist
key_ref <- paste(snp$bp, pmin(snp$a1, snp$a2), pmax(snp$a1, snp$a2))
key_hm3 <- paste(hm3$BP, pmin(hm3$A1, hm3$A2), pmax(hm3$A1, hm3$A2))
hm3_row <- match(key_ref, key_hm3)
keep <- which(!is.na(hm3_row))
keep <- keep[!duplicated(hm3_row[keep])]
if (length(keep) < 2) {
  saveRDS(empty, out_rds)
  quit(save = "no")
}

R <- res$cormat[keep, keep, drop = FALSE]
rows <- hm3_row[keep]
bp <- hm3$BP[rows]
ord <- order(bp)
R <- R[ord, ord, drop = FALSE]
rows <- rows[ord]
bp <- bp[ord]

ut <- which(upper.tri(R), arr.ind = TRUE)
a <- ut[, 1]
b <- ut[, 2]
dist <- bp[b] - bp[a]
in_anchor <- bp[a] >= win_start & bp[a] <= anchor_end
sel <- in_anchor & dist <= step_bp

out <- empty
out$pairs <- data.frame(i = rows[a[sel]], j = rows[b[sel]], r = R[ut[sel, , drop = FALSE]])
out$present <- rows[bp >= win_start & bp <= anchor_end]
saveRDS(out, out_rds)

message(sprintf("Kept %d of %d window SNPs as HapMap3; saved %d pairs in %s",
  length(rows), nrow(snp), nrow(out$pairs), format(Sys.time() - t0)))
