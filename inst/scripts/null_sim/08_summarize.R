#!/usr/bin/env Rscript
# Summarizes the null and non-null metrics (plan section 7): paired contrasts of
# PRS-MIX at the generating weight against each single-population arm, with
# the simulated GWAS replicate as the unit, plus weight-sweep curves.
#
# Each contrast is computed within replicate (same simulated GWAS and target
# individuals for both arms), as a log ratio for variances and errors and a
# plain difference for LD-fit scores, r^2 and sign agreement. Reported: mean
# over replicates, bootstrap 95% CI (resampling replicates), paired t-test p,
# and Holm-adjusted p across the primary contrasts within each phi setting.
# Primary contrasts (fixed in the plan before the full run): LD fit (main
# delta), null score variance and non-null r^2 in each target ancestry, each
# against both single-population arms.
#
# Usage:
#   Rscript 08_summarize.R <metrics_dir> [w0=0.7] [delta=0.001] [pops=EUR,ASN] [n_boot=10000]
#
# Writes <metrics_dir>/summary_contrasts.tsv and summary_curves.tsv.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args("Rscript 08_summarize.R <metrics_dir> [w0=0.7] [delta=0.001] [pops=EUR,ASN] [n_boot=10000]", 1)
mdir <- args[[1]]
w0 <- as.numeric(ns_arg(args, 2, "0.7"))
delta0 <- as.numeric(ns_arg(args, 3, "0.001"))
pops <- toupper(strsplit(ns_arg(args, 4, "EUR,ASN"), ",")[[1]])
n_boot <- as.integer(ns_arg(args, 5, "10000"))
singles <- paste0("single_", pops)
mix0 <- sprintf("mix_w%.1f", w0)

fit <- ns_read_tsv(file.path(mdir, "ld_fit.tsv"))
null <- ns_read_tsv(file.path(mdir, "prscs_metrics.tsv"))
nonnull <- if (file.exists(file.path(mdir, "nonnull_metrics.tsv"))) ns_read_tsv(file.path(mdir, "nonnull_metrics.tsv")) else NULL

set.seed(1)
contrast <- function(d, label) {
  d <- d[is.finite(d)]
  n <- length(d)
  boot <- replicate(n_boot, mean(sample(d, n, replace = TRUE)))
  data.frame(contrast = label, n_rep = n, mean = mean(d),
             ci_lo = stats::quantile(boot, 0.025, names = FALSE),
             ci_hi = stats::quantile(boot, 0.975, names = FALSE),
             p = if (n > 1 && stats::sd(d) > 0) stats::t.test(d)$p.value else NA_real_)
}

# Paired values of `col` for arm a vs arm b within one phi setting.
paired <- function(tab, phi_mode, col, a, b, log_ratio) {
  x <- tab[tab$phi_mode == phi_mode, ]
  va <- x[x$arm == a, c("rep", col)]
  vb <- x[x$arm == b, c("rep", col)]
  m <- merge(va, vb, by = "rep")
  if (log_ratio) log(m[[2]] / m[[3]]) else m[[2]] - m[[3]]
}

rows <- list()
add <- function(phi_mode, family, metric, scale, single, d) {
  r <- contrast(d, sprintf("%s: %s vs %s", metric, mix0, single))
  rows[[length(rows) + 1]] <<- cbind(data.frame(phi_mode = phi_mode, family = family, metric = metric,
                                                scale = scale, vs = single), r[, -1])
}

# LD fit does not involve PRS-CS, so it is the same for both phi settings.
for (dl in sort(unique(fit$delta))) {
  for (k in seq_along(pops)) {
    a <- fit[fit$w == w0 & fit$delta == dl, c("rep", "L")]
    b <- fit[fit$w == (if (k == 1) 1 else 0) & fit$delta == dl, c("rep", "L")]
    m <- merge(a, b, by = "rep")
    for (pm in unique(null$phi_mode)) {
      add(pm, if (dl == delta0) "primary" else "secondary", sprintf("ld_fit_L_delta%g", dl), "difference",
          singles[k], m$L.x - m$L.y)
    }
  }
}

for (pm in unique(null$phi_mode)) {
  for (s in singles) {
    for (p in pops) add(pm, "primary", paste0("null_score_var_", p), "log_ratio", s,
                        paired(null, pm, paste0("score_var_", p), mix0, s, TRUE))
    add(pm, "secondary", "null_sum_beta2", "log_ratio", s, paired(null, pm, "sum_beta2", mix0, s, TRUE))
    add(pm, "secondary", "null_abs_beta_q99", "log_ratio", s, paired(null, pm, "abs_beta_q99", mix0, s, TRUE))
    if (pm == "auto") add(pm, "secondary", "null_phi", "log_ratio", s, paired(null, pm, "phi", mix0, s, TRUE))
    if (!is.null(nonnull)) {
      for (p in pops) add(pm, "primary", paste0("nonnull_r2_", p), "difference", s,
                          paired(nonnull, pm, paste0("r2_", p), mix0, s, FALSE))
      for (p in pops) add(pm, "secondary", paste0("nonnull_mse_gv_", p), "log_ratio", s,
                          paired(nonnull, pm, paste0("mse_gv_", p), mix0, s, TRUE))
      add(pm, "secondary", "nonnull_err_rmix", "log_ratio", s, paired(nonnull, pm, "err_rmix", mix0, s, TRUE))
      add(pm, "secondary", "nonnull_err_plain", "log_ratio", s, paired(nonnull, pm, "err_plain", mix0, s, TRUE))
      add(pm, "secondary", "nonnull_sign_agree", "difference", s, paired(nonnull, pm, "sign_agree", mix0, s, FALSE))
      if (pm == "auto") add(pm, "secondary", "nonnull_phi", "log_ratio", s, paired(nonnull, pm, "phi", mix0, s, TRUE))
    }
  }
}
res <- do.call(rbind, rows)
res$p_holm <- NA_real_
for (pm in unique(res$phi_mode)) {
  i <- res$phi_mode == pm & res$family == "primary"
  res$p_holm[i] <- stats::p.adjust(res$p[i], "holm")
}
res$pct <- ifelse(res$scale == "log_ratio", 100 * (exp(res$mean) - 1), NA_real_)
ns_write_tsv(res, file.path(mdir, "summary_contrasts.tsv"))

# Weight-sweep curves: mean and bootstrap CI per weight, over the mix_w arms.
curve_rows <- list()
curve <- function(tab, phi_mode, col, source) {
  x <- tab[tab$phi_mode == phi_mode & grepl("^mix_w", tab$arm), ]
  if (!any(is.finite(x[[col]]))) return(invisible(NULL))  # e.g. phi under fixed-phi runs
  for (wv in sort(unique(x$w))) {
    r <- contrast(x[[col]][x$w == wv], col)
    curve_rows[[length(curve_rows) + 1]] <<- data.frame(source = source, phi_mode = phi_mode, metric = col,
                                                        w = wv, mean = r$mean, ci_lo = r$ci_lo, ci_hi = r$ci_hi)
  }
}
for (dl in sort(unique(fit$delta))) {
  x <- fit[fit$delta == dl, ]
  ref <- stats::setNames(x$L[x$w == w0], x$rep[x$w == w0])
  x$dL <- x$L - ref[as.character(x$rep)]
  for (wv in sort(unique(x$w))) {
    r <- contrast(x$dL[x$w == wv], "dL")
    curve_rows[[length(curve_rows) + 1]] <- data.frame(source = "null", phi_mode = "none",
                                                       metric = sprintf("ld_fit_dL_delta%g", dl), w = wv,
                                                       mean = r$mean, ci_lo = r$ci_lo, ci_hi = r$ci_hi)
  }
}
for (pm in unique(null$phi_mode)) {
  for (col in c(paste0("score_var_", pops), "sum_beta2", "phi")) curve(null, pm, col, "null")
  if (!is.null(nonnull)) for (col in c(paste0("r2_", pops), "err_rmix", "sign_agree", "phi")) curve(nonnull, pm, col, "nonnull")
}
curves <- do.call(rbind, curve_rows)
curves <- curves[is.finite(curves$mean), ]
ns_write_tsv(curves, file.path(mdir, "summary_curves.tsv"))
message(sprintf("Wrote %d contrasts and %d curve points to %s", nrow(res), nrow(curves), mdir))
