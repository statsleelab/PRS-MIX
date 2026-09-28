#!/usr/bin/env Rscript
# Simulates a null GWAS summary-statistics file: per-SNP z-scores with no true
# genotype-phenotype association, but with the same LD-induced correlation
# structure real finite-sample GWAS noise would have. Under the null,
# Cov(z_i, z_j) ~= r_ij (the reference-panel LD between i and j), not zero --
# drawing z ~ MVN(0, R_block) per LD block reproduces that, rather than the
# unrealistic (and much easier to get "right" by accident) i.i.d.-per-SNP
# alternative. This is the standard construction used to calibration-test
# PRS methods (PRS-CS, LDpred) against different LD reference panels.
#
# Each block's own shrunk, positive-definite ldblk matrix (as already written
# by write_prscs_ld_chr()/generate_prscs_ld_panels()) is treated as the "true"
# population LD for simulation purposes -- so the same panel_dir argument
# supplies both the correlation structure to simulate under and the allele
# frequencies used to convert the standardized z-score onto the per-allele
# BETA scale PRS-CS expects.
#
# Usage:
#   Rscript simulate_null_gwas.R <panel_dir> <chr> <n_gwas> <output_sumstats> [seed=1]
#
#   <panel_dir>        LD panel directory (e.g. ldblk_1kg_eur33kg) containing
#                      ldblk_1kg_chr<chr>.hdf5 and snpinfo_1kg_hm3.
#   <chr>              Chromosome to simulate.
#   <n_gwas>           GWAS sample size to simulate under (sets the
#                      per-SNP standard error via se = 1/sqrt(2*N*maf*(1-maf))).
#   <output_sumstats>  Where to write the null sumstats file (SNP/A1/A2/BETA/P,
#                      tab-delimited -- same format as a real PRS-CS sumstats file).
#   [seed]             RNG seed, defaults to 1. Vary this across replicates to
#                      assess effect-size stability under repeated null draws.
#
# Example:
#   Rscript inst/scripts/simulate_null_gwas.R \
#     gauss3_maf_test_40blk/maf_0p01/aligned/ldblk_1kg_eur33kg \
#     22 112311 null_sumstats_chr22_seed1.tsv 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript simulate_null_gwas.R <panel_dir> <chr> <n_gwas> <output_sumstats> [seed=1]")
}
panel_dir <- args[[1]]
chr <- as.integer(args[[2]])
n_gwas <- as.numeric(args[[3]])
output_path <- args[[4]]
seed <- if (length(args) >= 5) as.integer(args[[5]]) else 1L

if (!requireNamespace("hdf5r", quietly = TRUE)) {
  stop("Package 'hdf5r' is required. Please install it with install.packages('hdf5r').")
}

set.seed(seed)

snpinfo <- utils::read.table(file.path(panel_dir, "snpinfo_1kg_hm3"), header = TRUE, stringsAsFactors = FALSE)
maf_by_snp <- stats::setNames(snpinfo$MAF, snpinfo$SNP)

h5_path <- file.path(panel_dir, sprintf("ldblk_1kg_chr%d.hdf5", chr))
h5 <- hdf5r::H5File$new(h5_path, mode = "r")
on.exit(h5$close_all())

blk_names <- grep("^blk_", names(h5), value = TRUE)
blk_names <- blk_names[order(as.integer(sub("blk_", "", blk_names)))]

message(sprintf("Simulating null GWAS: chr%d, %d blocks, n_gwas=%s, seed=%d",
  chr, length(blk_names), format(n_gwas, big.mark = ","), seed))

out_parts <- vector("list", length(blk_names))

for (i in seq_along(blk_names)) {
  blk <- h5[[blk_names[i]]]
  m_snps <- blk$attr_open("m_snps")$read()
  if (m_snps < 2) next

  snplist <- blk[["snplist"]]$read()
  a1 <- blk[["a1"]]$read()
  a2 <- blk[["a2"]]$read()
  R <- as.matrix(blk[["ldblk"]]$read())

  maf <- maf_by_snp[snplist]
  if (anyNA(maf)) {
    stop(sprintf("Block %s: %d SNP(s) missing from snpinfo MAF lookup.", blk_names[i], sum(is.na(maf))))
  }

  # z ~ MVN(0, R): draw iid standard normal w, apply R's upper-triangular
  # Cholesky factor U (R = t(U) %*% U) so that Cov(t(U) %*% w) == R.
  U <- chol(R)
  w <- stats::rnorm(m_snps)
  z <- as.vector(crossprod(U, w))

  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  se <- 1 / sqrt(2 * n_gwas * maf * (1 - maf))
  beta <- z * se

  out_parts[[i]] <- data.frame(
    SNP = snplist, A1 = a1, A2 = a2, BETA = beta, P = p,
    stringsAsFactors = FALSE
  )

  if (i %% 20 == 0 || i == length(blk_names)) {
    message(sprintf("  block %d/%d done", i, length(blk_names)))
  }
}

out <- do.call(rbind, out_parts)
utils::write.table(out, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("Wrote %d null SNPs to %s", nrow(out), output_path))
