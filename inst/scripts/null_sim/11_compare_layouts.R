#!/usr/bin/env Rscript
# Compares two runs of the null/non-null simulation that differ only in the
# panel block layout (e.g. EUR-derived cut points vs 250 kb tiles). The runs
# must share the sample split, analysis SNPs and simulated GWAS (same seeds),
# so each replicate is the same GWAS under both layouts and differences are
# paired by replicate.
#
# For each phi setting and arm (EUR-only, PRS-MIX at w0, ASN-only), reports the
# mean of each metric under both layouts and the paired difference (b - a) with
# a normal-approximation 95% CI over replicates.
#
# Usage:
#   Rscript 11_compare_layouts.R <metrics_dir_a> <metrics_dir_b> <out_tsv> \
#     [label_a=eurcuts] [label_b=fixed250kb] [w0=0.7]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 11_compare_layouts.R <metrics_dir_a> <metrics_dir_b> <out_tsv> ",
       "[label_a=eurcuts] [label_b=fixed250kb] [w0=0.7]", call. = FALSE)
}
arg <- function(i, default) if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
dir_a <- args[[1]]
dir_b <- args[[2]]
out_tsv <- args[[3]]
label_a <- arg(4, "eurcuts")
label_b <- arg(5, "fixed250kb")
w0 <- as.numeric(arg(6, "0.7"))

arms <- c("single_EUR", sprintf("mix_w%.1f", w0), "single_ASN")
tables <- list(
  nonnull = c("r2_EUR", "r2_ASN", "mse_gv_EUR", "mse_gv_ASN", "err_rmix", "sign_agree"),
  prscs = c("score_var_EUR", "score_var_ASN", "sum_beta2")
)
files <- c(nonnull = "nonnull_metrics.tsv", prscs = "prscs_metrics.tsv")

rows <- list()
for (tab in names(tables)) {
  a <- utils::read.delim(file.path(dir_a, files[[tab]]))
  b <- utils::read.delim(file.path(dir_b, files[[tab]]))
  for (phi in unique(a$phi_mode)) for (arm in arms) for (m in tables[[tab]]) {
    x <- a[a$phi_mode == phi & a$arm == arm, c("rep", m)]
    y <- b[b$phi_mode == phi & b$arm == arm, c("rep", m)]
    z <- merge(x, y, by = "rep", suffixes = c(".a", ".b"))
    d <- z[[3]] - z[[2]]
    se <- stats::sd(d) / sqrt(length(d))
    rows[[length(rows) + 1]] <- data.frame(
      phi_mode = phi, arm = arm, metric = m, n_rep = length(d),
      mean_a = mean(z[[2]]), mean_b = mean(z[[3]]),
      diff = mean(d), ci_lo = mean(d) - 1.96 * se, ci_hi = mean(d) + 1.96 * se,
      pct = 100 * mean(d) / mean(z[[2]])
    )
  }
}
out <- do.call(rbind, rows)
names(out)[names(out) == "mean_a"] <- paste0("mean_", label_a)
names(out)[names(out) == "mean_b"] <- paste0("mean_", label_b)
utils::write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %d rows (%s vs %s) to %s", nrow(out), label_b, label_a, out_tsv))
