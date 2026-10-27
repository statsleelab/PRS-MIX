#!/usr/bin/env Rscript
# Null simulation, step 4: simulate null GWAS on the discovery genotypes and
# combine ancestries as a sample-size-weighted (Stouffer) meta-analysis.
#
# For replicate r and ancestry k (plan section 5.1):
#   y_k ~ N(0, I)                     a phenotype unrelated to genotype
#   z_k = Z_k' y_k / sqrt(n_k)        Cov(z_k) = Z_k'Z_k / n_k, the discovery LD
#   z   = sum_k sqrt(w_k) z_k         Cov(z)   = sum_k w_k Z_k'Z_k / n_k
# with w the generating composition (default 70:30 EUR:ASN). The genotypes are
# never block-diagonal, so LD across block boundaries is kept.
#
# Sumstats are written as BETA = z / sqrt(N), SE = 1 / sqrt(N), so PRS-CS's
# BETA/SE parser recovers exactly z / sqrt(N) at --n_gwas=N.
#
# Replicate r draws its phenotypes from seed + r, so replicates can be added
# later without changing earlier ones.
#
# Usage:
#   Rscript 04_simulate_null.R <geno_dir> <out_dir> <first_rep> <last_rep> \
#     [seed=1000] [w=EUR:0.7,ASN:0.3] [n_gwas=200000]
#
# Outputs in <out_dir>: sumstats/rep_<r>.tsv, z/rep_<r>.rds (named z vector),
# params.tsv.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args(paste(
  "Rscript 04_simulate_null.R <geno_dir> <out_dir> <first_rep> <last_rep>",
  "[seed=1000] [w=EUR:0.7,ASN:0.3] [n_gwas=200000]"
), 4)
geno_dir <- args[[1]]
out_dir <- args[[2]]
reps <- seq(as.integer(args[[3]]), as.integer(args[[4]]))
seed <- as.integer(ns_arg(args, 5, "1000"))
w_spec <- strsplit(strsplit(ns_arg(args, 6, "EUR:0.7,ASN:0.3"), ",")[[1]], ":")
w <- stats::setNames(as.numeric(vapply(w_spec, `[`, "", 2)), toupper(vapply(w_spec, `[`, "", 1)))
n_gwas <- as.numeric(ns_arg(args, 7, "200000"))
if (any(w <= 0) || abs(sum(w) - 1) > 1e-8) stop("weights must be positive and sum to 1", call. = FALSE)

snps <- ns_read_tsv(file.path(geno_dir, "analysis_snps.tsv"))
dir.create(file.path(out_dir, "sumstats"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "z"), recursive = TRUE, showWarnings = FALSE)

Z <- lapply(stats::setNames(names(w), names(w)), function(pop) {
  x <- ns_load_standardized(geno_dir, pop, "discovery")
  stopifnot(identical(colnames(x), snps$SNP))
  x
})
for (pop in names(Z)) message(sprintf("%s discovery: %d individuals x %d SNPs", pop, nrow(Z[[pop]]), ncol(Z[[pop]])))

ns_write_tsv(data.frame(param = c("seed", "n_gwas", paste0("w_", names(w)), paste0("n_discovery_", names(w))),
                        value = c(seed, n_gwas, w, vapply(Z, nrow, integer(1)))),
             file.path(out_dir, "params.tsv"))

for (r in reps) {
  set.seed(seed + r)
  z <- numeric(nrow(snps))
  for (pop in names(w)) {
    y <- stats::rnorm(nrow(Z[[pop]]))
    z <- z + sqrt(w[[pop]]) * as.vector(crossprod(Z[[pop]], y)) / sqrt(nrow(Z[[pop]]))
  }
  names(z) <- snps$SNP
  saveRDS(z, file.path(out_dir, "z", sprintf("rep_%03d.rds", r)))
  ns_write_tsv(data.frame(SNP = snps$SNP, A1 = snps$A1, A2 = snps$A2,
                          BETA = sprintf("%.10g", z / sqrt(n_gwas)), SE = sprintf("%.10g", 1 / sqrt(n_gwas))),
               file.path(out_dir, "sumstats", sprintf("rep_%03d.tsv", r)))
  message(sprintf("rep %d: mean z %.3f, var z %.3f, max |z| %.2f", r, mean(z), stats::var(z), max(abs(z))))
}
