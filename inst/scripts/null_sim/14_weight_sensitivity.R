#!/usr/bin/env Rscript
# Weight sensitivity across the ancestry-pair x GWAS-composition sweep, from the
# existing per-replicate metrics (no new PRS-CS runs).
#
# Scope: exploratory, not curated paper evidence. Chromosome 22 only, fixed
# 250 kb LD blocks, PRS-CS with automatic phi, default sampler settings, the
# first max_rep (default 20) replicates of each of the 12 settings (4 ancestry
# pairs x GWAS compositions 90:10, 70:30, 50:50). For each setting the
# generating weight w0 is the weight on the first-listed ancestry A; B is the
# second. PRS-MIX was run at every w in 0.0, 0.1, ..., 1.0 (w = 1 is A-only and
# w = 0 is B-only PRS-CS). Everything is paired within replicate.
#
#   1. Misspecification (non-null GWAS): r^2 with the true genetic value and
#      MSE at every weight, relative to PRS-MIX at w0. The relative r^2 change
#      is computed per replicate as r^2(w) / r^2(w0) - 1, so tol is a fraction
#      of the generating-weight r^2 (tol = 0.02: a loss of 2% of r^2(w0), e.g.
#      0.006 when r^2(w0) = 0.3), not an absolute r^2 difference. A weight is
#      called acceptable when the lower 95% bound of this change is at least
#      -tol; the threshold was fixed before the curves were examined. Intervals are pointwise t intervals over replicates and are
#      exploratory: nothing is adjusted for the many settings and weights.
#   2. Paired contrasts at w0: PRS-MIX minus A-only and minus B-only PRS-CS, in
#      r^2 (difference) and MSE (geometric-mean percentage change).
#   3. Null LD-fit recovery: the weight minimising the Gaussian LD-fit loss
#      (ridge 0.001) per replicate. The loss uses the NULL GWAS z-scores, which
#      are generated separately from the non-null GWAS, so this asks only
#      whether LD fit recovers the composition; it is not a weight selector
#      validated for prediction.
#   4. Oracle (retrospective): the weight with the highest mean target r^2 and
#      the better single panel, both chosen using the true genetic values of the
#      same targets. They are benchmarks a user could not reproduce.
#
#   weight_sensitivity.html  Tables 1-4 and two figures.
#   weight_sensitivity.tsv   Every number, unrounded, with n and 95% interval.
#   weight_curves_r2.pdf/.png, weight_curves_mse.pdf/.png
#
# Usage:
#   Rscript inst/scripts/null_sim/14_weight_sensitivity.R <root> <out_dir> [max_rep=20] [tol=0.02]
#   <root> holds <pair>/w<w>/metrics/{nonnull_metrics,ld_fit}.tsv.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages(library(ggplot2))
args <- ns_args("Rscript 14_weight_sensitivity.R <root> <out_dir> [max_rep=20] [tol=0.02]", 2)
root <- args[[1]]
out_dir <- args[[2]]
max_rep <- as.integer(ns_arg(args, 3, "20"))
tol <- as.numeric(ns_arg(args, 4, "0.02"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pair_order <- c("eur_asn", "eur_sas", "eur_afr", "asn_afr")
dirs <- unlist(lapply(file.path(root, pair_order), function(p) {
  if (!dir.exists(p)) return(character(0))
  w <- list.files(p, pattern = "^w[0-9.]+$", full.names = TRUE, include.dirs = TRUE)
  w[order(-as.numeric(sub("^w", "", basename(w))))]
}))
weights <- round(seq(1, 0, by = -0.1), 1)
arm_of <- function(w) sprintf("mix_w%.1f", w)

per_rep <- function(tab, arm, col) {
  x <- tab[tab$arm == arm, ]
  stats::setNames(x[[col]], x$rep)
}
paired <- function(a, b, f = `-`) {
  k <- intersect(names(a), names(b))
  f(a[k], b[k])
}
cells <- list()
put <- function(analysis, s, target, w, column, x, scale = "mean") {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
  q <- if (n > 1) stats::qt(0.975, n - 1) else NA_real_
  lo <- m - q * se
  hi <- m + q * se
  if (scale == "log_pct" && n) {
    lo <- 100 * (exp(lo) - 1)
    hi <- 100 * (exp(hi) - 1)
    se <- 100 * exp(m) * se
    m <- 100 * (exp(m) - 1)
  }
  cells[[length(cells) + 1]] <<- data.frame(analysis = analysis, pair = s$pair, mix = s$mix, w0 = s$w0, target = target,
                                            w = w, column = column, n = n, mean = m, se = se, ci_lo = lo, ci_hi = hi)
}

settings <- list()
for (d in dirs) {
  pops <- toupper(strsplit(basename(dirname(d)), "_")[[1]])
  A <- pops[1]
  B <- pops[2]
  w0 <- as.numeric(sub("^w", "", basename(d)))
  m <- file.path(d, "metrics")
  nn <- ns_read_tsv(file.path(m, "nonnull_metrics.tsv"))
  ld <- ns_read_tsv(file.path(m, "ld_fit.tsv"))
  nn <- nn[nn$phi_mode == "auto" & nn$rep <= max_rep, ]
  ld <- ld[abs(ld$delta - 0.001) < 1e-12 & ld$rep <= max_rep, ]
  stopifnot(all(arm_of(weights) %in% nn$arm), all(weights %in% round(ld$w, 1)))
  s <- list(pair = sprintf("%s + %s", A, B), mix = sprintf("%.0f:%.0f", 100 * w0, 100 * (1 - w0)), A = A, B = B, w0 = w0)
  settings[[length(settings) + 1]] <- s
  ref <- arm_of(w0)

  for (tg in c(A, B)) {
    r2 <- paste0("r2_", tg)
    mse <- paste0("mse_gv_", tg)
    # 1. every weight against the matching weight
    for (w in weights) {
      put("curve", s, tg, w, "r2", per_rep(nn, arm_of(w), r2))
      put("curve", s, tg, w, "r2_rel", 100 * paired(per_rep(nn, arm_of(w), r2), per_rep(nn, ref, r2), function(a, b) a / b - 1))
      put("curve", s, tg, w, "mse_rel", paired(log(per_rep(nn, arm_of(w), mse)), log(per_rep(nn, ref, mse))), "log_pct")
    }
    # 2. PRS-MIX at w0 against each single panel
    for (sp in c(1, 0)) {
      who <- if (sp == 1) "A" else "B"
      put("contrast", s, tg, w0, paste0("d_r2_", who), paired(per_rep(nn, ref, r2), per_rep(nn, arm_of(sp), r2)))
      put("contrast", s, tg, w0, paste0("mse_", who), paired(log(per_rep(nn, ref, mse)), log(per_rep(nn, arm_of(sp), mse))), "log_pct")
    }
  }
  # 3. null LD-fit argmin per replicate
  amin <- vapply(split(ld, ld$rep), function(x) round(x$w[which.min(x$L)], 1), numeric(1))
  put("ld_recovery", s, "", NA, "argmin_w", amin)
  put("ld_recovery", s, "", NA, "abs_err", abs(amin - w0))
  put("ld_recovery", s, "", NA, "exact", 100 * (abs(amin - w0) < 1e-9))
  put("ld_recovery", s, "", NA, "within_0.1", 100 * (abs(amin - w0) < 0.1 + 1e-9))
}
cells <- do.call(rbind, cells)

# 4. oracle, from the mean curves (retrospective)
oracle <- list()
for (s in settings) for (tg in c(s$A, s$B)) {
  cv <- cells[cells$analysis == "curve" & cells$pair == s$pair & cells$mix == s$mix & cells$target == tg & cells$column == "r2", ]
  best <- cv[which.max(cv$mean), ]
  single <- cv[cv$w %in% c(1, 0), ]
  bs <- single[which.max(single$mean), ]
  oracle[[length(oracle) + 1]] <- data.frame(
    analysis = "oracle", pair = s$pair, mix = s$mix, w0 = s$w0, target = tg,
    w = c(s$w0, best$w, bs$w), column = c("r2_matching", "r2_best_weight", "r2_best_single"),
    n = c(cv$n[cv$w == s$w0], best$n, bs$n), mean = c(cv$mean[cv$w == s$w0], best$mean, bs$mean),
    se = c(cv$se[cv$w == s$w0], best$se, bs$se), ci_lo = NA_real_, ci_hi = NA_real_)
}
cells <- rbind(cells, do.call(rbind, oracle))
ns_write_tsv(cells, file.path(out_dir, "weight_sensitivity.tsv"))

# ---- Figures ---------------------------------------------------------------
ink <- "#1b1b1b"
ink2 <- "#5b5f66"
col <- c("#1f5a96", "#c4552d")
pair_levels <- unique(vapply(settings, `[[`, "", "pair"))
curve_fig <- function(column, ylab, band) {
  d <- cells[cells$analysis == "curve" & cells$column == column, ]
  d$role <- factor(ifelse(d$target == sub(" \\+.*", "", d$pair), "first", "second"), c("first", "second"),
                   c("Targets of the first-listed ancestry", "Targets of the second-listed ancestry"))
  d$pair <- factor(d$pair, pair_levels)
  d$mix <- factor(paste("GWAS", d$mix), paste("GWAS", c("90:10", "70:30", "50:50")))
  v <- unique(d[, c("pair", "mix", "w0")])
  ggplot(d, aes(w, mean, colour = role, fill = role)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2], fill = "grey88", alpha = 0.6) +
    geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
    geom_vline(data = v, aes(xintercept = w0), linetype = "22", colour = ink2, linewidth = 0.3) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.9) +
    facet_grid(pair ~ mix, scales = "free_y") +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_colour_manual(values = col, name = NULL) +
    scale_fill_manual(values = col, name = NULL) +
    labs(x = "PRS-MIX weight on the first-listed ancestry (0 = second-ancestry PRS-CS, 1 = first-ancestry PRS-CS)", y = ylab) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid = element_blank(), axis.line = element_line(colour = ink2, linewidth = 0.3),
          axis.ticks = element_line(colour = ink2, linewidth = 0.3), axis.text = element_text(colour = ink2),
          axis.title = element_text(colour = ink), strip.text = element_text(face = "bold", colour = ink, size = 8),
          legend.position = "bottom", panel.spacing = grid::unit(9, "pt"))
}
figs <- list(
  r2 = curve_fig("r2_rel", "r² relative to PRS-MIX at the generating weight (%)", c(-100 * tol, 100 * tol)),
  mse = curve_fig("mse_rel", "MSE relative to PRS-MIX at the generating weight (%)", c(-100 * tol, 100 * tol)))
for (k in names(figs)) {
  ggsave(file.path(out_dir, sprintf("weight_curves_%s.pdf", k)), figs[[k]], width = 18, height = 19, units = "cm", device = grDevices::pdf)
  ggsave(file.path(out_dir, sprintf("weight_curves_%s.png", k)), figs[[k]], width = 18, height = 19, units = "cm", dpi = 300, bg = "white")
}

# ---- HTML ------------------------------------------------------------------
num <- function(x, d, signed = FALSE) {
  if (is.na(x)) return("&ndash;")
  s <- formatC(abs(x), format = "f", digits = d)
  if (as.numeric(s) == 0) return(s)
  paste0(if (x < 0) "&minus;" else if (signed) "+" else "", s)
}
get <- function(analysis, s, tg, column, w = NULL) {
  r <- cells[cells$analysis == analysis & cells$pair == s$pair & cells$mix == s$mix & cells$target == tg & cells$column == column, ]
  if (!is.null(w)) r <- r[!is.na(r$w) & abs(r$w - w) < 1e-9, ]
  if (nrow(r)) r[1, ] else NULL
}
ci_cell <- function(r, d, signed = TRUE, flag = FALSE) {
  if (is.null(r)) return("&ndash;")
  out <- sprintf("%s <span class=\"ci\">[%s, %s]</span>", num(r$mean, d, signed), num(r$ci_lo, d, signed), num(r$ci_hi, d, signed))
  if (flag && r$ci_lo < -100 * tol) out <- paste0("<b>", out, "</b>")
  out
}
html_table <- function(groups, heads, row_fn, per_target = TRUE) {
  top <- paste0("<th></th><th></th>", if (per_target) "<th></th>" else "",
                paste(sprintf("<th colspan=\"%d\" class=\"grp\">%s</th>", vapply(groups, `[[`, 0, 2), vapply(groups, `[[`, "", 1)), collapse = ""))
  sub <- paste0("<th>Pair</th><th>GWAS mix</th>", if (per_target) "<th>Targets</th>" else "", paste(sprintf("<th>%s</th>", heads), collapse = ""))
  body <- character(0)
  last_pair <- ""
  for (s in settings) {
    tgs <- if (per_target) c(s$A, s$B) else ""
    for (i in seq_along(tgs)) {
      first_of_pair <- s$pair != last_pair && i == 1
      cls <- if (first_of_pair && last_pair != "") " class=\"newpair\"" else ""
      body <- c(body, sprintf("<tr%s><td>%s</td><td>%s</td>%s%s</tr>", cls, if (first_of_pair) s$pair else "",
                              if (i == 1) s$mix else "", if (per_target) sprintf("<td>%s</td>", tgs[i]) else "",
                              paste(sprintf("<td>%s</td>", row_fn(s, tgs[i])), collapse = "")))
      last_pair <- s$pair
    }
  }
  c("<table>", "<thead>", sprintf("<tr class=\"top\">%s</tr>", top), sprintf("<tr class=\"sub\">%s</tr>", sub), "</thead>",
    "<tbody>", body, "</tbody>", "</table>")
}
offsets <- c(-0.2, -0.1, 0.1, 0.2)
t_mis <- html_table(list(list("Relative <i>r</i><sup>2</sup> change (%) when the weight is off by", 4)),
                    sprintf("%+.1f", offsets),
                    function(s, tg) vapply(offsets, function(o) {
                      w <- round(s$w0 + o, 1)
                      if (w < 0 || w > 1) "&ndash;" else ci_cell(get("curve", s, tg, "r2_rel", w), 1, flag = TRUE)
                    }, ""))
t_con <- html_table(list(list("PRS-MIX minus first-ancestry PRS-CS", 2), list("PRS-MIX minus second-ancestry PRS-CS", 2)),
                    c("<i>r</i><sup>2</sup>", "MSE (%)", "<i>r</i><sup>2</sup>", "MSE (%)"),
                    function(s, tg) c(ci_cell(get("contrast", s, tg, "d_r2_A"), 4), ci_cell(get("contrast", s, tg, "mse_A"), 1),
                                      ci_cell(get("contrast", s, tg, "d_r2_B"), 4), ci_cell(get("contrast", s, tg, "mse_B"), 1)))
t_ld <- html_table(list(list("Weight minimising null LD-fit loss, per replicate", 4)),
                   c("Mean (SE)", "Mean absolute error", "Exactly w<sub>0</sub> (%)", "Within &plusmn;0.1 (%)"),
                   function(s, tg) {
                     a <- get("ld_recovery", s, "", "argmin_w")
                     c(sprintf("%s (%s)", num(a$mean, 2), num(a$se, 2)), num(get("ld_recovery", s, "", "abs_err")$mean, 2),
                       num(get("ld_recovery", s, "", "exact")$mean, 0), num(get("ld_recovery", s, "", "within_0.1")$mean, 0))
                   }, per_target = FALSE)
t_or <- html_table(list(list("PRS-MIX at w<sub>0</sub>", 1), list("Best weight (oracle)", 2), list("Better single panel (oracle)", 2)),
                   c("<i>r</i><sup>2</sup>", "Weight", "<i>r</i><sup>2</sup>", "Panel", "<i>r</i><sup>2</sup>"),
                   function(s, tg) {
                     a <- get("oracle", s, tg, "r2_matching")
                     b <- get("oracle", s, tg, "r2_best_weight")
                     c2 <- get("oracle", s, tg, "r2_best_single")
                     c(sprintf("%s (%s)", num(a$mean, 3), num(a$se, 3)), sprintf("%.1f", b$w), sprintf("%s (%s)", num(b$mean, 3), num(b$se, 3)),
                       if (c2$w == 1) s$A else s$B, sprintf("%s (%s)", num(c2$mean, 3), num(c2$se, 3)))
                   })
common <- paste("Exploratory analysis, not curated paper evidence. Chromosome 22, fixed 250 kb LD blocks, PRS-CS with automatic &phi; and default sampler settings; replicates 1&ndash;", max_rep,
                " per setting. w<sub>0</sub> is the generating weight on the first-listed ancestry; PRS-MIX at w = 1 and 0 is",
                " first- and second-ancestry PRS-CS. Square brackets are pointwise 95% t intervals over replicates, paired within",
                " replicate; they are exploratory and not adjusted for the number of settings and weights.", sep = "")
sec <- function(title, caption, body) c("<section>", sprintf("<h2>%s</h2>", title), sprintf("<p class=\"cap\">%s</p>", caption), body, "</section>")
wrap <- function(tab) c("<div class=\"wrap\">", tab, "</div>")
img <- function(k, alt) sprintf("<figure><img alt=\"%s\" src=\"data:image/png;base64,%s\"></figure>", alt,
                                base64enc::base64encode(file.path(out_dir, sprintf("weight_curves_%s.png", k))))
html <- c(
  "<!doctype html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">", "<title>Weight sensitivity</title>",
  "<!-- Generated by inst/scripts/null_sim/14_weight_sensitivity.R; do not edit. -->", "<style>",
  "body { margin: 0; padding: 40px 24px 56px; background: #fff; color: #111; font: 15px/1.45 Georgia, \"Times New Roman\", Times, serif; }",
  "main { max-width: 1180px; margin: 0 auto; display: flex; flex-direction: column; gap: 44px; }",
  "h2 { font: 600 17px/1.3 -apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif; margin: 0 0 8px; }",
  ".cap { margin: 0 0 14px; font-size: 14px; }", ".wrap { overflow-x: auto; }", ".ci { color: #555; font-size: 12.5px; }",
  "table { border-collapse: collapse; width: 100%; font-variant-numeric: lining-nums tabular-nums; font-size: 14px; }",
  "th, td { padding: 3px 10px; text-align: center; white-space: nowrap; vertical-align: baseline; }",
  "th:nth-child(-n+3), td:nth-child(-n+3) { text-align: left; }",
  "thead th { font-weight: 700; }", "thead tr.top th { border-top: 2px solid #111; padding-top: 6px; }",
  "thead tr.sub th { border-bottom: 1px solid #111; padding-bottom: 6px; }",
  "thead th.grp { border-bottom: 1px solid #111; }", "tbody tr.newpair td { border-top: 1px solid #111; }",
  "tbody tr:last-child td { border-bottom: 2px solid #111; padding-bottom: 6px; }",
  "figure { margin: 0; } figure img { width: 100%; max-width: 900px; display: block; }",
  "footer { font: 12px/1.5 -apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif; color: #666; }",
  "</style>", "</head>", "<body>", "<main>",
  sec("Table 1. Weight misspecification (non-null GWAS)",
      paste0("Per-replicate percentage change in <i>r</i><sup>2</sup> with the true genetic value when PRS-MIX uses a weight off from w<sub>0</sub> by the column offset, relative to PRS-MIX at w<sub>0</sub>; negative values mean lower accuracy. The change is relative, <i>r</i><sup>2</sup>(w)/<i>r</i><sup>2</sup>(w<sub>0</sub>) &minus; 1 per replicate, so the threshold is a fraction of the generating-weight <i>r</i><sup>2</sup>, not an absolute difference. Bold: the lower 95% bound falls below &minus;", 100 * tol,
             "%, the acceptable relative loss fixed before the curves were examined; unbolded cells are within it. &ndash;: the weight is outside [0, 1]. ", common),
      wrap(t_mis)),
  sec("Table 2. PRS-MIX at the generating weight against the single panels (non-null GWAS)",
      paste("Paired differences in <i>r</i><sup>2</sup> (positive favours PRS-MIX) and geometric-mean percentage change in mean squared error against the true genetic value (negative favours PRS-MIX).", common),
      wrap(t_con)),
  sec("Table 3. Recovery of the composition by null LD fit",
      paste("Per replicate, the weight in 0, 0.1, &hellip;, 1 that minimises the Gaussian LD-fit loss (ridge 0.001) of the mixed panel against the null GWAS z-scores; mean (Monte Carlo SE) over replicates. The null GWAS is simulated separately from the non-null GWAS in Tables 1, 2 and 4, so this asks only whether LD fit recovers the composition; it does not validate LD fit as a way to choose the weight for prediction.", common),
      wrap(t_ld)),
  sec("Table 4. Oracle benchmarks (retrospective)",
      paste("The weight with the highest mean <i>r</i><sup>2</sup>, and the better of the two single panels, both chosen with the true genetic values of the same target individuals; mean (Monte Carlo SE). These choices are not available to a user and are optimistic by construction; they show how much a perfectly informed choice could gain over the generating weight.", common),
      wrap(t_or)),
  sec("Figure 1. <i>r</i><sup>2</sup> across weights",
      paste0("Mean percentage change in <i>r</i><sup>2</sup> relative to PRS-MIX at the generating weight (dashed vertical line), with pointwise 95% intervals; grey band &plusmn;", 100 * tol, "%. Free y-axis per row."),
      img("r2", "r squared relative to the generating weight across PRS-MIX weights, by setting")),
  sec("Figure 2. Mean squared error across weights",
      paste0("Geometric-mean percentage change in MSE against the true genetic value relative to PRS-MIX at the generating weight, with pointwise 95% intervals; negative is better; grey band &plusmn;", 100 * tol, "%. Free y-axis per row."),
      img("mse", "MSE relative to the generating weight across PRS-MIX weights, by setting")),
  sprintf("<footer>Generated from <code>%s</code> by <code>inst/scripts/null_sim/14_weight_sensitivity.R</code> (replicates 1&ndash;%d, tol = %g) at PRS-MIX %s. Every number is in <code>weight_sensitivity.tsv</code>.</footer>",
          root, max_rep, tol, substr(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), 1, 7)),
  "</main>", "</body>", "</html>")
writeLines(html, file.path(out_dir, "weight_sensitivity.html"))
message(sprintf("Wrote %s: %d settings, %d cells, replicates 1-%d", out_dir, length(settings), nrow(cells), max_rep))
