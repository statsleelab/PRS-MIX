#!/usr/bin/env Rscript
# Non-null simulation metrics (plan section 7.2), per replicate, phi setting
# and arm, from the true effects written by 05_simulate_nonnull.R and each
# arm's PRS-CS output (standardized betas, --beta_std=TRUE).
#
#   err_rmix     ||bhat - beta||^2_R / ||beta||^2_R with R = sum_k w_k Z_k'Z_k / n_k,
#                the generating covariance -- the same yardstick for every arm
#   err_plain    ||bhat - beta||^2 / ||beta||^2
#   sign_agree   fraction of causal SNPs with sign(bhat) = sign(beta)
#   per target ancestry k (Z_T standardized with discovery constants; g = Z_T beta,
#   s = Z_T bhat):
#     r2_<POP>       cor(s, g)^2, recovery of genetic value (not phenotype R^2)
#     mse_gv_<POP>   mean of ((s - mean s) - (g - mean g))^2, scale as well as order
#     vs single_EUR and single_ASN: cor_<REF>_<POP>, varratio_<REF>_<POP>,
#     top10_<REF>_<POP> (overlap of the top score decile), top05_ and top20_
#     (top 5% and 20%), spear_<REF>_<POP> (Spearman rank correlation)
#   phi            posterior global shrinkage estimate (auto runs only)
#
# Usage:
#   Rscript 07_metrics_nonnull.R <sim_dir> <chr> [pops=EUR,ASN]
#
# Reads <sim_dir>/{geno,nonnull,prscs_nonnull}; writes
# <sim_dir>/metrics/{nonnull_metrics,nonnull_identity}.tsv.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args("Rscript 07_metrics_nonnull.R <sim_dir> <chr> [pops=EUR,ASN]", 2)
sim_dir <- args[[1]]
chr <- as.integer(args[[2]])
pops <- toupper(strsplit(ns_arg(args, 3, "EUR,ASN"), ",")[[1]])
stopifnot(length(pops) == 2)

geno_dir <- file.path(sim_dir, "geno")
out_dir <- file.path(sim_dir, "metrics")
dir.create(out_dir, showWarnings = FALSE)
snps <- ns_read_tsv(file.path(geno_dir, "analysis_snps.tsv"))

params <- ns_read_tsv(file.path(sim_dir, "nonnull", "params.tsv"))
w <- stats::setNames(params$value[match(paste0("w_", pops), params$param)], pops)
stopifnot(!anyNA(w))

disc <- lapply(stats::setNames(pops, pops), function(p) ns_load_standardized(geno_dir, p, "discovery"))
targets <- lapply(stats::setNames(pops, pops), function(p) ns_load_standardized(geno_dir, p, "target"))

rmix_norm2 <- function(b) sum(vapply(pops, function(p) w[[p]] * sum(as.vector(disc[[p]] %*% b)^2) / nrow(disc[[p]]), numeric(1)))
top_overlap <- function(a, b, frac = 0.1) {
  k <- ceiling(round(length(a) * frac, 6))
  length(intersect(order(-a)[seq_len(k)], order(-b)[seq_len(k)])) / k
}

refs <- paste0("single_", pops)
met <- list()
ident <- list()
for (phi_dir in list.dirs(file.path(sim_dir, "prscs_nonnull"), recursive = FALSE)) {
  phi_mode <- sub("^phi_", "", basename(phi_dir))
  for (rep_dir in list.dirs(phi_dir, recursive = FALSE)) {
    rep <- as.integer(sub("^rep_", "", basename(rep_dir)))
    beta <- readRDS(file.path(sim_dir, "nonnull", "beta", sprintf("rep_%03d.rds", rep)))
    stopifnot(identical(names(beta), snps$SNP))
    causal <- beta != 0
    beta_norm2 <- rmix_norm2(beta)
    gv <- lapply(targets, function(x) as.vector(x %*% beta))

    files <- list.files(rep_dir, sprintf("_pst_eff_.*_chr%d\\.txt$", chr), full.names = TRUE)
    arms <- sub("_pst_eff_.*$", "", basename(files))
    betas <- stats::setNames(lapply(files, ns_read_beta, snps = snps), arms)
    scores <- lapply(betas, function(b) lapply(targets, function(x) as.vector(x %*% b)))

    for (arm in arms) {
      bh <- betas[[arm]]
      row <- data.frame(phi_mode = phi_mode, rep = rep, arm = arm, w = ns_arm_weight(arm, pops),
                        phi = ns_read_phi(file.path(rep_dir, paste0(arm, ".log"))),
                        err_rmix = rmix_norm2(bh - beta) / beta_norm2,
                        err_plain = sum((bh - beta)^2) / sum(beta^2),
                        sign_agree = mean(sign(bh[causal]) == sign(beta[causal])))
      for (p in pops) {
        s <- scores[[arm]][[p]]
        g <- gv[[p]]
        row[[paste0("r2_", p)]] <- if (stats::sd(s) > 0) stats::cor(s, g)^2 else NA_real_
        row[[paste0("mse_gv_", p)]] <- mean(((s - mean(s)) - (g - mean(g)))^2)
        for (ref in refs) {
          tag <- paste0(sub("^single_", "", ref), "_", p)
          r_s <- if (ref %in% arms) scores[[ref]][[p]] else NULL
          row[[paste0("cor_", tag)]] <- if (is.null(r_s)) NA_real_ else stats::cor(s, r_s)
          row[[paste0("varratio_", tag)]] <- if (is.null(r_s)) NA_real_ else stats::var(s) / stats::var(r_s)
          row[[paste0("top10_", tag)]] <- if (is.null(r_s)) NA_real_ else top_overlap(s, r_s)
          row[[paste0("top05_", tag)]] <- if (is.null(r_s)) NA_real_ else top_overlap(s, r_s, 0.05)
          row[[paste0("top20_", tag)]] <- if (is.null(r_s)) NA_real_ else top_overlap(s, r_s, 0.2)
          row[[paste0("spear_", tag)]] <- if (is.null(r_s)) NA_real_ else stats::cor(s, r_s, method = "spearman")
        }
      }
      met[[length(met) + 1]] <- row
    }
    for (pr in list(c("mix_w1.0", refs[1]), c("mix_w0.0", refs[2]))) {
      if (all(pr %in% arms)) {
        ident[[length(ident) + 1]] <- data.frame(phi_mode = phi_mode, rep = rep, mix = pr[1], single = pr[2],
                                                 max_abs_diff = max(abs(betas[[pr[1]]] - betas[[pr[2]]])))
      }
    }
  }
}
met <- do.call(rbind, met)
ns_write_tsv(met, file.path(out_dir, "nonnull_metrics.tsv"))
ident <- if (length(ident)) do.call(rbind, ident) else data.frame()
ns_write_tsv(ident, file.path(out_dir, "nonnull_identity.tsv"))
message(sprintf("Wrote %d non-null metric rows; identity check max |diff| = %s",
                nrow(met), if (nrow(ident)) format(max(ident$max_abs_diff)) else "NA"))
