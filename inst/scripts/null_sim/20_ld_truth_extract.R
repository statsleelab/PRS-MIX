#!/usr/bin/env Rscript
# Null simulation: extract the LD needed to compare each candidate LD matrix
# with the LD the null GWAS was actually simulated under. Runs where the
# discovery genotypes are (RedHawk); 21_ld_truth_figures.R draws the figures
# from its outputs.
#
# For ancestries A and B and GWAS composition w (weight on A), 04_simulate_null.R
# gives Cov(z) = w R_A,disc + (1 - w) R_B,disc, with R_k,disc = Z_k'Z_k / n_k
# from the discovery individuals (the "true GWAS LD"). The candidates are built
# from the separate reference individuals' panels:
#   PRS-MIX   w R_A,panel + (1 - w) R_B,panel
#   A alone   R_A,panel  (standard PRS-CS with the majority-ancestry panel)
#   B alone   R_B,panel
# Neither R_k,disc nor the panels depend on w, so one extraction serves every
# composition of the pair (all compositions share one split, panels and genotypes).
#
# Outputs in <out_dir>:
#   <pair>_block_scan.tsv   every block with >= 2 analysis SNPs, for w = 0.9,
#                           0.7, 0.5: mean |true - candidate| and RMS difference
#                           over the block's off-diagonal SNP pairs. Used to
#                           choose and justify the drawn block.
#   <pair>_block<b>_ld.tsv  for each requested block, every SNP pair (i <= j):
#                           discovery LD and panel LD of A and B.
#
# Usage (from the PRS-MIX root, or with the path to this script):
#   Rscript inst/scripts/null_sim/20_ld_truth_extract.R <sim_dir> <out_dir> <pops e.g. EUR,ASN> <blocks e.g. 75> [layout=fixed250kb] [chr=22]
#
#   <sim_dir>  a null-simulation directory with panels/<layout>/ and geno/
#              (e.g. .../outputs/sim250/eur_asn/w0.7).

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
args <- ns_args("Rscript 20_ld_truth_extract.R <sim_dir> <out_dir> <pops> <blocks> [layout=fixed250kb] [chr=22]", 4)
sim <- args[[1]]
out_dir <- args[[2]]
pops <- toupper(strsplit(args[[3]], ",")[[1]])
blocks <- as.integer(strsplit(args[[4]], ",")[[1]])
layout <- ns_arg(args, 5, "fixed250kb")
chr <- as.integer(ns_arg(args, 6, "22"))
stopifnot(length(pops) == 2)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pair <- tolower(paste(pops, collapse = "_"))
ws <- c(0.9, 0.7, 0.5)

an <- ns_read_tsv(file.path(sim, "geno", "analysis_snps.tsv"))
pa <- ns_read_panel(ns_panel_dir(file.path(sim, "panels", layout), pops[1]), chr)
pb <- ns_read_panel(ns_panel_dir(file.path(sim, "panels", layout), pops[2]), chr)
stopifnot(length(pa) == length(pb))
X <- lapply(stats::setNames(pops, pops), function(p) {
  x <- ns_load_standardized(file.path(sim, "geno"), p, "discovery")
  stopifnot(identical(colnames(x), an$SNP))
  x
})

# Panel block i restricted to the analysis SNPs, with the discovery LD of both ancestries.
block_ld <- function(i) {
  stopifnot(identical(pa[[i]]$snps, pb[[i]]$snps))
  k <- which(pa[[i]]$snps %in% an$SNP)
  if (length(k) < 2) return(NULL)
  s <- pa[[i]]$snps[k]
  j <- match(s, an$SNP)
  disc <- lapply(X, function(x) crossprod(x[, j, drop = FALSE]) / nrow(x))
  list(snps = s, bp = an$BP[j], da = disc[[1]], db = disc[[2]],
       ra = pa[[i]]$ld[k, k, drop = FALSE], rb = pb[[i]]$ld[k, k, drop = FALSE])
}

scan <- list()
for (i in seq_along(pa)) {
  b <- block_ld(i)
  if (is.null(b)) next
  ut <- upper.tri(b$ra)
  for (w in ws) {
    tr <- w * b$da + (1 - w) * b$db
    e <- list(prsmix = (tr - (w * b$ra + (1 - w) * b$rb))[ut], a_only = (tr - b$ra)[ut], b_only = (tr - b$rb)[ut])
    scan[[length(scan) + 1]] <- data.frame(
      pair = paste(pops, collapse = "+"), w = w, block = i, n_snps = length(b$snps),
      bp_start = min(b$bp), bp_end = max(b$bp), mean_abs_panel_diff = mean(abs((b$ra - b$rb)[ut])),
      mad_prsmix = mean(abs(e$prsmix)), mad_a_only = mean(abs(e$a_only)), mad_b_only = mean(abs(e$b_only)),
      rms_prsmix = sqrt(mean(e$prsmix^2)), rms_a_only = sqrt(mean(e$a_only^2)), rms_b_only = sqrt(mean(e$b_only^2)))
  }
  if (i %in% blocks) {
    n <- length(b$snps)
    keep <- upper.tri(b$ra, diag = TRUE)
    ii <- row(b$ra)[keep]; jj <- col(b$ra)[keep]
    ns_write_tsv(data.frame(i = ii, j = jj, snp_i = b$snps[ii], snp_j = b$snps[jj], bp_i = b$bp[ii], bp_j = b$bp[jj],
                            r_disc_a = b$da[keep], r_disc_b = b$db[keep], r_panel_a = b$ra[keep], r_panel_b = b$rb[keep]),
                 file.path(out_dir, sprintf("%s_block%d_ld.tsv", pair, i)))
  }
}
scan <- do.call(rbind, scan)
ns_write_tsv(scan, file.path(out_dir, sprintf("%s_block_scan.tsv", pair)))
ns_write_tsv(data.frame(param = c("sim_dir", "pops", "layout", "chr", "blocks", sprintf("n_discovery_%s", pops)),
                        value = c(sim, paste(pops, collapse = ","), layout, chr, paste(blocks, collapse = ","),
                                  vapply(X, nrow, integer(1)))),
             file.path(out_dir, sprintf("%s_params.tsv", pair)))
message(sprintf("%s: %d blocks scanned; wrote LD for block(s) %s", paste(pops, collapse = "+"),
                length(unique(scan$block)), paste(blocks, collapse = ",")))
