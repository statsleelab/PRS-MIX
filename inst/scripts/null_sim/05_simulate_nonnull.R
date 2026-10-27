#!/usr/bin/env Rscript
# Non-null simulation (plan sections 3.3 and 5.2): the null construction of
# 04_simulate_null.R plus a planted signal, at the same generating composition.
#
# For replicate r, with standardized discovery genotypes Z_k (n_k x m) and
# composition weights w_k:
#   beta     sparse: a random p_causal fraction of SNPs get N(0, 1) effects,
#            rescaled so that sum_k w_k ||Z_k beta||^2 / n_k = v_g, i.e.
#            beta' R_mix beta = v_g under the generating covariance
#   z        sqrt(N) sum_k w_k R_k beta + sum_k sqrt(w_k) Z_k' y_k / sqrt(n_k),
#            R_k beta = Z_k'(Z_k beta) / n_k,  y_k ~ N(0, I)
# This is z ~ N(sqrt(N) R_mix beta, R_mix): a sample-size-weighted
# meta-analysis of ancestry GWAS with nominal sizes w_k N and shared
# standardized effects. Each replicate draws new causal SNPs, effects and
# noise from seed + r.
#
# v_g is a signal variance in units of residual variance one, not an exact
# phenotypic heritability (plan section 5.2).
#
# Usage:
#   Rscript 05_simulate_nonnull.R <geno_dir> <out_dir> <first_rep> <last_rep> \
#     [seed=5000] [w=EUR:0.7,ASN:0.3] [n_gwas=200000] [p_causal=0.01] [v_g=0.01]
#
# Outputs in <out_dir>: sumstats/rep_<r>.tsv, z/rep_<r>.rds, beta/rep_<r>.rds
# (named true standardized effects), params.tsv.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args(paste(
  "Rscript 05_simulate_nonnull.R <geno_dir> <out_dir> <first_rep> <last_rep>",
  "[seed=5000] [w=EUR:0.7,ASN:0.3] [n_gwas=200000] [p_causal=0.01] [v_g=0.01]"
), 4)
geno_dir <- args[[1]]
out_dir <- args[[2]]
reps <- seq(as.integer(args[[3]]), as.integer(args[[4]]))
seed <- as.integer(ns_arg(args, 5, "5000"))
w_spec <- strsplit(strsplit(ns_arg(args, 6, "EUR:0.7,ASN:0.3"), ",")[[1]], ":")
w <- stats::setNames(as.numeric(vapply(w_spec, `[`, "", 2)), toupper(vapply(w_spec, `[`, "", 1)))
n_gwas <- as.numeric(ns_arg(args, 7, "200000"))
p_causal <- as.numeric(ns_arg(args, 8, "0.01"))
v_g <- as.numeric(ns_arg(args, 9, "0.01"))
if (any(w <= 0) || abs(sum(w) - 1) > 1e-8) stop("weights must be positive and sum to 1", call. = FALSE)
if (!(p_causal > 0 && p_causal <= 1) || !(v_g > 0)) stop("need 0 < p_causal <= 1 and v_g > 0", call. = FALSE)

snps <- ns_read_tsv(file.path(geno_dir, "analysis_snps.tsv"))
m <- nrow(snps)
n_causal <- max(1L, as.integer(round(p_causal * m)))
for (d in c("sumstats", "z", "beta")) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)

Z <- lapply(stats::setNames(names(w), names(w)), function(pop) {
  x <- ns_load_standardized(geno_dir, pop, "discovery")
  stopifnot(identical(colnames(x), snps$SNP))
  x
})

ns_write_tsv(data.frame(
  param = c("seed", "n_gwas", "p_causal", "n_causal", "v_g", paste0("w_", names(w)), paste0("n_discovery_", names(w))),
  value = c(seed, n_gwas, p_causal, n_causal, v_g, w, vapply(Z, nrow, integer(1)))
), file.path(out_dir, "params.tsv"))

for (r in reps) {
  set.seed(seed + r)
  beta <- numeric(m)
  causal <- sort(sample.int(m, n_causal))
  beta[causal] <- stats::rnorm(n_causal)

  # Genetic values in each ancestry's discovery set, used for both scaling and signal.
  gv <- lapply(Z, function(x) as.vector(x %*% beta))
  scale <- sqrt(v_g / sum(vapply(names(w), function(p) w[[p]] * sum(gv[[p]]^2) / nrow(Z[[p]]), numeric(1))))
  beta <- beta * scale

  signal <- numeric(m)
  noise <- numeric(m)
  for (pop in names(w)) {
    n <- nrow(Z[[pop]])
    signal <- signal + w[[pop]] * as.vector(crossprod(Z[[pop]], gv[[pop]] * scale)) / n
    noise <- noise + sqrt(w[[pop]]) * as.vector(crossprod(Z[[pop]], stats::rnorm(n))) / sqrt(n)
  }
  z <- sqrt(n_gwas) * signal + noise
  names(z) <- names(beta) <- snps$SNP

  saveRDS(z, file.path(out_dir, "z", sprintf("rep_%03d.rds", r)))
  saveRDS(beta, file.path(out_dir, "beta", sprintf("rep_%03d.rds", r)))
  ns_write_tsv(data.frame(SNP = snps$SNP, A1 = snps$A1, A2 = snps$A2,
                          BETA = sprintf("%.10g", z / sqrt(n_gwas)), SE = sprintf("%.10g", 1 / sqrt(n_gwas))),
               file.path(out_dir, "sumstats", sprintf("rep_%03d.tsv", r)))
  message(sprintf("rep %d: %d causal SNPs, mean chi2 %.3f, max |z| %.2f, max |z| at causal %.2f",
                  r, n_causal, mean(z^2), max(abs(z)), max(abs(z[causal]))))
}
