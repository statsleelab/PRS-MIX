#!/usr/bin/env Rscript
# Internal results page (tables only) for the null and non-null simulations,
# from the metric tables of 06/07 and the contrasts of 08_summarize.R.
#
#   1. Run settings and checks (replicates, identity of PRS-MIX endpoints with
#      vanilla PRS-CS).
#   2. Each metric for EUR-only, PRS-MIX at the generating weight and ASN-only
#      (mean over replicates, 95% bootstrap CI).
#   3. Paired contrasts of PRS-MIX against each single-population arm.
#   4. Every metric across the weight sweep (means).
#
# Usage:
#   Rscript 10_report.R <metrics_dir> [out_html=<metrics_dir>/results_summary.html] [w0=0.7] [delta=0.001]

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args("Rscript 10_report.R <metrics_dir> [out_html] [w0=0.7] [delta=0.001]", 1)
mdir <- args[[1]]
out_html <- ns_arg(args, 2, file.path(mdir, "results_summary.html"))
w0 <- as.numeric(ns_arg(args, 3, "0.7"))
delta0 <- as.numeric(ns_arg(args, 4, "0.001"))
pops <- c("EUR", "ASN")
mix0 <- sprintf("mix_w%.1f", w0)

fit <- ns_read_tsv(file.path(mdir, "ld_fit.tsv"))
null <- ns_read_tsv(file.path(mdir, "prscs_metrics.tsv"))
nonnull <- ns_read_tsv(file.path(mdir, "nonnull_metrics.tsv"))
con <- ns_read_tsv(file.path(mdir, "summary_contrasts.tsv"))
ident <- rbind(cbind(kind = "null", ns_read_tsv(file.path(mdir, "identity.tsv"))),
               cbind(kind = "non-null", ns_read_tsv(file.path(mdir, "nonnull_identity.tsv"))))

esc <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE)
html_table <- function(df, cls = "") {
  head <- paste0("<tr>", paste0("<th>", esc(names(df)), "</th>", collapse = ""), "</tr>")
  body <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table class=\"", cls, "\">\n", head, "\n", paste(body, collapse = "\n"), "\n</table>")
}
set.seed(1)
boot_ci <- function(x, n_boot = 10000) {
  x <- x[is.finite(x)]
  b <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
  c(mean(x), stats::quantile(b, c(0.025, 0.975), names = FALSE))
}
fmt_ci <- function(v, digits) {
  if (any(!is.finite(v))) return("&ndash;")
  sprintf("%s <span class=\"ci\">[%s, %s]</span>",
          formatC(v[1], format = "f", digits = digits), formatC(v[2], format = "f", digits = digits),
          formatC(v[3], format = "f", digits = digits))
}

# --- 2. Arm values ------------------------------------------------------------
arm_defs <- c("EUR-only" = "single_EUR", "PRS-MIX (w = 0.7)" = mix0, "ASN-only" = "single_ASN")
arm_defs <- stats::setNames(arm_defs, sub("0.7", format(w0), names(arm_defs), fixed = TRUE))
metric_defs <- list(
  list("LD-fit score (δ = %g), lower = better", "fit", "L", 1, 3),
  list("Null: score variance, EUR patients (×10⁻⁴), lower = better", "null", "score_var_EUR", 1e4, 3),
  list("Null: score variance, ASN patients (×10⁻⁴), lower = better", "null", "score_var_ASN", 1e4, 3),
  list("Null: estimated φ (×10⁻⁴)", "null", "phi", 1e4, 3),
  list("Non-null: r² with genetic value, EUR patients, higher = better", "nonnull", "r2_EUR", 1, 3),
  list("Non-null: r² with genetic value, ASN patients, higher = better", "nonnull", "r2_ASN", 1, 3),
  list("Non-null: effect error (R_mix-weighted), lower = better", "nonnull", "err_rmix", 1, 3),
  list("Non-null: sign agreement at causal SNPs, higher = better", "nonnull", "sign_agree", 1, 3),
  list("Non-null: estimated φ (×10⁻⁴)", "nonnull", "phi", 1e4, 3))

arm_value <- function(src, col, arm, phi_mode) {
  if (src == "fit") {
    w <- if (arm == "single_EUR") 1 else if (arm == "single_ASN") 0 else w0
    return(fit$L[fit$delta == delta0 & abs(fit$w - w) < 1e-9])
  }
  tab <- if (src == "null") null else nonnull
  tab[[col]][tab$phi_mode == phi_mode & tab$arm == arm]
}
arm_table <- function(phi_mode) {
  rows <- lapply(metric_defs, function(m) {
    if (m[[3]] == "phi" && phi_mode != "auto") return(NULL)
    vals <- vapply(arm_defs, function(a) fmt_ci(boot_ci(arm_value(m[[2]], m[[3]], a, phi_mode)) * m[[4]], m[[5]]), "")
    c(Metric = if (m[[2]] == "fit") sprintf(m[[1]], delta0) else m[[1]], vals)
  })
  as.data.frame(do.call(rbind, rows), check.names = FALSE)
}

# --- 3. Contrasts -------------------------------------------------------------
contrast_table <- function(phi_mode) {
  x <- con[con$phi_mode == phi_mode, ]
  x <- x[order(x$family != "primary", x$metric, x$vs), ]
  est <- function(m, lo, hi, scale) {
    if (scale == "log_ratio") sprintf("%+.1f%% <span class=\"ci\">[%+.1f, %+.1f]</span>",
                                      100 * (exp(m) - 1), 100 * (exp(lo) - 1), 100 * (exp(hi) - 1))
    else sprintf("%+.4f <span class=\"ci\">[%+.4f, %+.4f]</span>", m, lo, hi)
  }
  p_fmt <- function(p) ifelse(is.na(p), "&ndash;", formatC(p, format = "g", digits = 2))
  labels <- c(
    null_score_var_EUR = "Null score variance, EUR patients", null_score_var_ASN = "Null score variance, ASN patients",
    nonnull_r2_EUR = "Non-null r² with genetic value, EUR patients", nonnull_r2_ASN = "Non-null r² with genetic value, ASN patients",
    null_sum_beta2 = "Null sum of β²", null_abs_beta_q99 = "Null 99th percentile of |β|",
    null_phi = "Null estimated φ", nonnull_phi = "Non-null estimated φ",
    nonnull_err_rmix = "Non-null effect error (R_mix-weighted)", nonnull_err_plain = "Non-null effect error (plain)",
    nonnull_mse_gv_EUR = "Non-null genetic-value MSE, EUR patients", nonnull_mse_gv_ASN = "Non-null genetic-value MSE, ASN patients",
    nonnull_sign_agree = "Non-null sign agreement at causal SNPs")
  lab <- ifelse(grepl("^ld_fit_L_delta", x$metric),
                sprintf("LD-fit score (δ = %s)", sub("^ld_fit_L_delta", "", x$metric)), labels[x$metric])
  data.frame(
    Family = x$family,
    Metric = ifelse(is.na(lab), x$metric, lab),
    `Compared with` = sub("^single_", "", paste0(x$vs, "-only")),
    `PRS-MIX minus single panel` = mapply(est, x$mean, x$ci_lo, x$ci_hi, x$scale),
    p = p_fmt(x$p),
    `Holm p` = p_fmt(x$p_holm),
    check.names = FALSE)
}

# --- 4. Sweep -----------------------------------------------------------------
sweep_table <- function(phi_mode) {
  ws <- round(seq(1, 0, by = -0.1), 1)
  row <- function(label, vals, digits) c(Metric = label, stats::setNames(formatC(vals, format = "f", digits = digits), paste0("w=", ws)))
  ld <- vapply(ws, function(w) mean(fit$L[fit$delta == delta0 & abs(fit$w - w) < 1e-9]), numeric(1))
  m <- function(tab, col) vapply(ws, function(w) mean(tab[[col]][tab$phi_mode == phi_mode & grepl("^mix_w", tab$arm) & abs(tab$w - w) < 1e-9]), numeric(1))
  rows <- list(
    row(sprintf("LD-fit score (δ = %g)", delta0), ld, 3),
    row("Null score var, EUR (×10⁻⁴)", 1e4 * m(null, "score_var_EUR"), 3),
    row("Null score var, ASN (×10⁻⁴)", 1e4 * m(null, "score_var_ASN"), 3),
    row("Non-null r², EUR", m(nonnull, "r2_EUR"), 3),
    row("Non-null r², ASN", m(nonnull, "r2_ASN"), 3),
    row("Non-null effect error", m(nonnull, "err_rmix"), 3))
  as.data.frame(do.call(rbind, rows), check.names = FALSE)
}

n_null <- length(unique(null$rep)); n_nonnull <- length(unique(nonnull$rep))
phi_modes <- c("auto", setdiff(unique(null$phi_mode), "auto"))
phi_label <- function(p) if (p == "auto") "PRS-CS auto φ" else paste0("PRS-CS φ = ", p)

checks <- data.frame(
  Check = c("Null replicates", "Non-null replicates", "PRS-CS runs (null + non-null)",
            "Max |β| difference, PRS-MIX endpoints vs vanilla PRS-CS"),
  Result = c(n_null, n_nonnull, nrow(null) + nrow(nonnull),
             sprintf("%g (%d comparisons)", max(ident$max_abs_diff), nrow(ident))))

sections <- c(
  "<h2>1. Run and checks</h2>",
  sprintf(paste0("<p>chr22, HapMap3 SNPs; panels built from the 50%% reference split of 33KG with EUR-derived block cut points; ",
                 "GWAS simulated on the 25%% discovery split as a %g:%g EUR:ASN sample-size-weighted meta-analysis at N = 200,000; ",
                 "scores evaluated in the 25%% target split. Non-null: 1%% causal SNPs, v<sub>g</sub> = 0.01. ",
                 "Each metric is computed per replicate and arm; CIs are 95%% bootstrap intervals over replicates.</p>"),
          100 * w0, 100 * (1 - w0)),
  html_table(checks))
for (pm in phi_modes) {
  sections <- c(sections,
    sprintf("<h2>2. Metric values by arm &mdash; %s</h2>", phi_label(pm)),
    "<p>The LD-fit score does not involve PRS-CS, so it is the same under every φ setting.</p>",
    html_table(arm_table(pm)),
    sprintf("<h2>3. Paired contrasts &mdash; %s</h2>", phi_label(pm)),
    paste0("<p>Within-replicate differences of PRS-MIX at w = ", w0, " minus the single-population arm: % change for variances ",
           "and errors (log ratio), plain difference otherwise. Holm correction across the primary contrasts within this φ setting. ",
           "With δ = 0 the single-population panels' near-singular blocks make their LD-fit scores very large, ",
           "so those differences reflect conditioning more than fit.</p>"),
    html_table(contrast_table(pm), "contrasts"),
    sprintf("<h2>4. Weight sweep (means) &mdash; %s</h2>", phi_label(pm)),
    "<p>w = EUR weight given to PRS-MIX; w = 1 and w = 0 equal EUR-only and ASN-only.</p>",
    html_table(sweep_table(pm), "sweep"))
}

page <- c(
  "<!DOCTYPE html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
  sprintf("<title>Null Simulation Results</title>"),
  "<style>",
  "body { font-family: -apple-system, BlinkMacSystemFont, \"Segoe UI\", sans-serif; max-width: 1100px; margin: 2rem auto; padding: 0 1rem; line-height: 1.5; color: #222; background: #fff; }",
  "h1 { border-bottom: 1px solid #ddd; padding-bottom: .3rem; } h2 { margin-top: 2rem; font-size: 1.15rem; }",
  "table { border-collapse: collapse; margin: .8rem 0; display: block; overflow-x: auto; font-size: .9rem; }",
  "th, td { border: 1px solid #ddd; padding: .3rem .6rem; text-align: left; vertical-align: top; white-space: nowrap; }",
  "th { background: #f6f8fa; } .ci { color: #777; font-size: .85em; }",
  "table.sweep td, table.sweep th { text-align: right; } table.sweep td:first-child, table.sweep th:first-child { text-align: left; }",
  "</style>", "</head>", "<body>",
  sprintf("<h1>Null and non-null simulation results (%d replicates)</h1>", max(n_null, n_nonnull)),
  "<p style=\"color:#666\">Internal summary, generated by inst/scripts/null_sim/10_report.R.</p>",
  sections, "</body>", "</html>")
writeLines(page, out_html, useBytes = FALSE)
message(sprintf("Wrote %s", out_html))
