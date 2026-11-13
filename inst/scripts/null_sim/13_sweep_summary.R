#!/usr/bin/env Rscript
# Summary tables and violin figure for the ancestry-pair x GWAS-composition
# sweep (one SIM directory per setting, <root>/<a>_<b>/w<w>/metrics/).
#
# For each setting the generating weight w0 is the weight on the first-listed
# ancestry A; B is the second. Methods: A-only PRS-CS (reference), its seed
# rerun, PRS-MIX at w0, and B-only PRS-CS. Every quantity is computed per
# replicate and paired within replicate (same simulated GWAS and target
# individuals), then summarised as mean (Monte Carlo SE) over replicates.
# Score-variance ratios use the log scale (geometric mean percentage change).
#
#   sweep_tables.html   Table 1 null GWAS, Table 2 non-null GWAS, Table 3 PRS-MIX
#                       against the seed rerun and the single panels.
#   sweep_tables.tsv    Every number in the tables, unrounded, with n.
#   sweep_violins.pdf/.png  Per-replicate distributions: (a) PRS-MIX minus seed
#                       rerun in top-10% overlap with A-only PRS-CS; (b) PRS-MIX
#                       minus A-only PRS-CS in r^2 with the true genetic value.
#
# Usage:
#   Rscript inst/scripts/null_sim/13_sweep_summary.R <root> <out_dir> [max_rep=20]
#   <root> holds <pair>/w<w>/metrics/{prscs_metrics,nonnull_metrics,ld_fit}.tsv.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages(library(ggplot2))
args <- ns_args("Rscript 13_sweep_summary.R <root> <out_dir> [max_rep=20]", 2)
root <- args[[1]]
out_dir <- args[[2]]
max_rep <- as.integer(ns_arg(args, 3, "20"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pair_order <- c("eur_asn", "eur_sas", "eur_afr", "asn_afr")
dirs <- list.files(root, pattern = "^w[0-9.]+$", full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
dirs <- unlist(lapply(file.path(root, pair_order), function(p) {
  if (!dir.exists(p)) return(character(0))
  w <- list.files(p, pattern = "^w[0-9.]+$", full.names = TRUE, include.dirs = TRUE)
  w[order(-as.numeric(sub("^w", "", basename(w))))]
}))

# ---- Per-replicate quantities ---------------------------------------------
per_rep <- function(tab, arm, col) {
  x <- tab[tab$arm == arm, ]
  stats::setNames(x[[col]], x$rep)
}
paired <- function(a, b) {
  k <- intersect(names(a), names(b))
  a[k] - b[k]
}
reps_long <- list()   # per-replicate values for the figure
cells <- list()
put <- function(table, s, target, column, x, scale = "mean") {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
  if (scale == "log_pct" && n) {
    se <- 100 * exp(m) * se
    m <- 100 * (exp(m) - 1)
  }
  cells[[length(cells) + 1]] <<- data.frame(table = table, pair = s$pair, mix = s$mix, target = target,
                                            column = column, n = n, mean = m, se = se)
}
settings <- list()
for (d in dirs) {
  pops <- toupper(strsplit(basename(dirname(d)), "_")[[1]])
  A <- pops[1]
  B <- pops[2]
  w0 <- as.numeric(sub("^w", "", basename(d)))
  m <- file.path(d, "metrics")
  pm <- ns_read_tsv(file.path(m, "prscs_metrics.tsv"))
  nn <- ns_read_tsv(file.path(m, "nonnull_metrics.tsv"))
  ld <- ns_read_tsv(file.path(m, "ld_fit.tsv"))
  pm <- pm[pm$phi_mode == "auto" & pm$rep <= max_rep, ]
  nn <- nn[nn$phi_mode == "auto" & nn$rep <= max_rep, ]
  ld <- ld[abs(ld$delta - 0.001) < 1e-12 & ld$rep <= max_rep, ]
  s <- list(pair = sprintf("%s + %s", A, B), mix = sprintf("%.0f:%.0f", 100 * w0, 100 * (1 - w0)), A = A, B = B, w0 = w0)
  settings[[length(settings) + 1]] <- s
  ref <- paste0("single_", A)
  seed <- paste0(ref, "_seed2")
  mix <- sprintf("mix_w%.1f", w0)
  other <- paste0("single_", B)
  ld_ex <- function(w) {
    a <- ld[abs(ld$w - w) < 1e-9, ]
    b <- ld[abs(ld$w - w0) < 1e-9, ]
    paired(stats::setNames(a$L, a$rep), stats::setNames(b$L, b$rep))
  }
  put("null", s, "", "ld_A", ld_ex(1))
  put("null", s, "", "ld_B", ld_ex(0))
  for (tg in c(A, B)) {
    lv <- function(arm, base) paired(log(per_rep(pm, arm, paste0("score_var_", tg))), log(per_rep(pm, base, paste0("score_var_", tg))))
    put("null", s, tg, "var_seed", lv(seed, ref), "log_pct")
    put("null", s, tg, "var_mix", lv(mix, ref), "log_pct")
    put("null", s, tg, "var_B", lv(other, ref), "log_pct")
    top <- sprintf("top10_%s_%s", A, tg)
    cor <- sprintf("cor_%s_%s", A, tg)
    r2 <- paste0("r2_", tg)
    put("nonnull", s, tg, "top_seed", 100 * per_rep(nn, seed, top))
    put("nonnull", s, tg, "top_mix", 100 * per_rep(nn, mix, top))
    put("nonnull", s, tg, "top_B", 100 * per_rep(nn, other, top))
    put("nonnull", s, tg, "r2_A", per_rep(nn, ref, r2))
    put("nonnull", s, tg, "r2_mix", per_rep(nn, mix, r2))
    put("nonnull", s, tg, "r2_B", per_rep(nn, other, r2))
    d_top <- 100 * paired(per_rep(nn, mix, top), per_rep(nn, seed, top))
    d_r2A <- paired(per_rep(nn, mix, r2), per_rep(nn, ref, r2))
    put("contrast", s, tg, "d_top_seed", d_top)
    put("contrast", s, tg, "d_cor_seed", 1000 * paired(per_rep(nn, mix, cor), per_rep(nn, seed, cor)))
    put("contrast", s, tg, "var_vs_seed", lv(mix, seed), "log_pct")
    put("contrast", s, tg, "d_r2_A", d_r2A)
    put("contrast", s, tg, "d_r2_B", paired(per_rep(nn, mix, r2), per_rep(nn, other, r2)))
    reps_long[[length(reps_long) + 1]] <- data.frame(
      pair = s$pair, mix = s$mix, target = tg, role = if (tg == A) "first" else "second",
      metric = rep(c("d_top", "d_r2"), c(length(d_top), length(d_r2A))), value = c(d_top, d_r2A))
  }
}
cells <- do.call(rbind, cells)
ns_write_tsv(cells, file.path(out_dir, "sweep_tables.tsv"))
reps_long <- do.call(rbind, reps_long)

# ---- Violin figure ---------------------------------------------------------
pair_levels <- vapply(settings, `[[`, "", "pair")
pair_levels <- unique(pair_levels)
reps_long$pair <- factor(reps_long$pair, pair_levels)
reps_long$mix <- factor(reps_long$mix, c("90:10", "70:30", "50:50"))
reps_long$role <- factor(reps_long$role, c("first", "second"),
                         c("Targets of the first-listed ancestry", "Targets of the second-listed ancestry"))
ink <- "#1b1b1b"
ink2 <- "#5b5f66"
fill <- c("#1f5a96", "#c4552d")
panel <- function(metric, ylab, title) {
  d <- reps_long[reps_long$metric == metric, ]
  ggplot(d, aes(x = mix, y = value, fill = role, colour = role)) +
    geom_hline(yintercept = 0, linetype = "22", colour = ink2, linewidth = 0.3) +
    geom_violin(position = position_dodge(0.8), width = 0.75, alpha = 0.25, linewidth = 0.3, scale = "width") +
    geom_boxplot(position = position_dodge(0.8), width = 0.18, outlier.shape = NA, alpha = 0.9,
                 fill = "white", linewidth = 0.35) +
    geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8, seed = 1),
               size = 0.55, alpha = 0.55, stroke = 0) +
    facet_wrap(~pair, nrow = 1) +
    scale_fill_manual(values = fill, name = NULL) +
    scale_colour_manual(values = fill, name = NULL) +
    labs(x = "GWAS composition (first:second ancestry)", y = ylab, title = title) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(panel.grid = element_blank(), axis.line = element_line(colour = ink2, linewidth = 0.3),
          axis.ticks = element_line(colour = ink2, linewidth = 0.3), axis.text = element_text(colour = ink2),
          axis.title = element_text(colour = ink), strip.text = element_text(face = "bold", colour = ink, size = 8),
          plot.title = element_text(face = "bold", size = 9, colour = ink), legend.position = "bottom",
          panel.spacing.x = grid::unit(10, "pt"))
}
pa <- panel("d_top", "PRS-MIX minus seed rerun\n(top-10% overlap, pp)", "a")
pb <- panel("d_r2", "PRS-MIX minus first-ancestry PRS-CS\n(r² with true genetic value)", "b")
fig <- patchwork::wrap_plots(pa, pb, ncol = 1, guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(out_dir, "sweep_violins.pdf"), fig, width = 18, height = 13, units = "cm", device = grDevices::pdf)
ggsave(file.path(out_dir, "sweep_violins.png"), fig, width = 18, height = 13, units = "cm", dpi = 300, bg = "white")

# ---- HTML tables -----------------------------------------------------------
num <- function(x, d, signed = FALSE) {
  s <- formatC(abs(x), format = "f", digits = d)
  if (as.numeric(s) == 0) return(s)
  paste0(if (x < 0) "&minus;" else if (signed) "+" else "", s)
}
cell <- function(table, s, tg, column, d, signed = FALSE) {
  r <- cells[cells$table == table & cells$pair == s$pair & cells$mix == s$mix & cells$target == tg & cells$column == column, ]
  if (nrow(r) == 0 || is.na(r$mean)) return("&ndash;")
  if (r$mean == 0 && isTRUE(r$se == 0)) return("0")
  unit <- 10^-d
  se <- if (is.na(r$se)) "" else if (r$se < unit / 2) sprintf(" (&lt;%s)", formatC(unit, format = "f", digits = d))
        else sprintf(" (%s)", formatC(r$se, format = "f", digits = d))
  paste0(num(r$mean, d, signed), se)
}
nrep <- function(s) {
  r <- cells[cells$pair == s$pair & cells$mix == s$mix & cells$column == "r2_A", "n"]
  r[1]
}
html_table <- function(groups, heads, row_fn, per_target = TRUE) {
  top <- paste0("<th></th><th></th><th></th>", if (per_target) "<th></th>" else "",
                paste(sprintf("<th colspan=\"%d\" class=\"grp\">%s</th>", vapply(groups, `[[`, 0, 2), vapply(groups, `[[`, "", 1)), collapse = ""))
  sub <- paste0("<th>Pair</th><th>GWAS mix</th>", if (per_target) "<th>Targets</th>" else "", "<th>Replicates</th>",
                paste(sprintf("<th>%s</th>", heads), collapse = ""))
  body <- character(0)
  last_pair <- ""
  for (s in settings) {
    tgs <- if (per_target) c(s$A, s$B) else ""
    for (i in seq_along(tgs)) {
      first_of_pair <- s$pair != last_pair && i == 1
      cls <- if (first_of_pair && last_pair != "") " class=\"newpair\"" else ""
      body <- c(body, sprintf("<tr%s><td>%s</td><td>%s</td>%s<td>%s</td>%s</tr>", cls,
                              if (first_of_pair) s$pair else "", if (i == 1) s$mix else "",
                              if (per_target) sprintf("<td>%s</td>", tgs[i]) else "",
                              if (i == 1) nrep(s) else "", paste(sprintf("<td>%s</td>", row_fn(s, tgs[i])), collapse = "")))
      last_pair <- s$pair
    }
  }
  c("<table>", "<thead>", sprintf("<tr class=\"top\">%s</tr>", top), sprintf("<tr class=\"sub\">%s</tr>", sub), "</thead>",
    "<tbody>", body, "</tbody>", "</table>")
}
t_null_ld <- html_table(list(list("LD-fit loss above the generating weight", 2)),
                        c("First-ancestry panel", "Second-ancestry panel"),
                        function(s, tg) c(cell("null", s, "", "ld_A", 3), cell("null", s, "", "ld_B", 3)), per_target = FALSE)
t_null <- html_table(list(list("Null score variance vs. first-ancestry PRS-CS (%)", 3)),
                     c("Seed rerun", "PRS-MIX", "Second-ancestry PRS-CS"),
                     function(s, tg) c(cell("null", s, tg, "var_seed", 1, TRUE), cell("null", s, tg, "var_mix", 1, TRUE),
                                       cell("null", s, tg, "var_B", 1, TRUE)))
t_nonnull <- html_table(list(list("Top-10% overlap with first-ancestry PRS-CS (%)", 3), list("<i>r</i><sup>2</sup> with true genetic value", 3)),
                        c("Seed rerun", "PRS-MIX", "Second-ancestry PRS-CS", "First-ancestry PRS-CS", "PRS-MIX", "Second-ancestry PRS-CS"),
                        function(s, tg) c(cell("nonnull", s, tg, "top_seed", 1), cell("nonnull", s, tg, "top_mix", 1),
                                          cell("nonnull", s, tg, "top_B", 1), cell("nonnull", s, tg, "r2_A", 3),
                                          cell("nonnull", s, tg, "r2_mix", 3), cell("nonnull", s, tg, "r2_B", 3)))
t_contrast <- html_table(list(list("PRS-MIX vs. seed rerun", 3), list("PRS-MIX minus single-panel PRS-CS (<i>r</i><sup>2</sup>)", 2)),
                         c("Top-10% overlap (pp)", "Correlation &times;10<sup>3</sup>", "Null score variance (%)",
                           "First ancestry", "Second ancestry"),
                         function(s, tg) c(cell("contrast", s, tg, "d_top_seed", 2, TRUE), cell("contrast", s, tg, "d_cor_seed", 2, TRUE),
                                           cell("contrast", s, tg, "var_vs_seed", 1, TRUE), cell("contrast", s, tg, "d_r2_A", 4, TRUE),
                                           cell("contrast", s, tg, "d_r2_B", 4, TRUE)))
common <- paste("Chromosome 22, fixed 250 kb LD blocks. Each setting is a GWAS combining two ancestries' discovery samples at the",
                "stated composition; PRS-MIX uses the matching weight. PRS-CS with automatic &phi;. Means over replicates with Monte",
                "Carlo standard errors in parentheses; within a replicate every method analyses the same simulated GWAS and scores",
                "the same target individuals. The seed rerun repeats first-ancestry PRS-CS with a different sampler seed and measures",
                "run-to-run variation.")
sec <- function(title, caption, tab) c("<section>", sprintf("<h2>%s</h2>", title), sprintf("<p class=\"cap\">%s</p>", caption),
                                       "<div class=\"wrap\">", tab, "</div>", "</section>")
fig_b64 <- base64enc::base64encode(file.path(out_dir, "sweep_violins.png"))
html <- c(
  "<!doctype html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">", "<title>sim250 sweep tables</title>",
  "<!-- Generated by inst/scripts/null_sim/13_sweep_summary.R; do not edit. -->", "<style>",
  "body { margin: 0; padding: 40px 24px 56px; background: #fff; color: #111; font: 15px/1.45 Georgia, \"Times New Roman\", Times, serif; }",
  "main { max-width: 1180px; margin: 0 auto; display: flex; flex-direction: column; gap: 44px; }",
  "h2 { font: 600 17px/1.3 -apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif; margin: 0 0 8px; }",
  ".cap { margin: 0 0 14px; font-size: 14px; }", ".wrap { overflow-x: auto; }",
  "table { border-collapse: collapse; width: 100%; font-variant-numeric: lining-nums tabular-nums; font-size: 14px; }",
  "th, td { padding: 3px 10px; text-align: center; white-space: nowrap; vertical-align: baseline; }",
  "th:nth-child(-n+3), td:nth-child(-n+3) { text-align: left; }",
  "thead th { font-weight: 700; }", "thead tr.top th { border-top: 2px solid #111; padding-top: 6px; }",
  "thead tr.sub th { border-bottom: 1px solid #111; padding-bottom: 6px; }",
  "thead th.grp { border-bottom: 1px solid #111; }", "tbody tr.newpair td { border-top: 1px solid #111; }",
  "tbody tr:last-child td { border-bottom: 2px solid #111; padding-bottom: 6px; }",
  "figure { margin: 0; } figure img { width: 100%; max-width: 1000px; display: block; }",
  "footer { font: 12px/1.5 -apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif; color: #666; }",
  "</style>", "</head>", "<body>", "<main>",
  sec("Table 1. Null GWAS simulation: LD fit", paste("Excess Gaussian LD-fit loss (ridge 0.001) of each single panel over the mixed panel at the generating weight, per SNP; lower is better, and the mixed panel is 0 by definition.", common), t_null_ld),
  sec("Table 2. Null GWAS simulation: score variance", paste("Variance of target-sample scores in a GWAS with no true effects, relative to first-ancestry PRS-CS for the same GWAS, as a geometric-mean percentage change; negative values mean less spread in scores that carry no true signal.", common), t_null),
  sec("Table 3. Non-null GWAS simulation", paste("Sparse causal effects shared across ancestries. Top-10% overlap is the share of target individuals in a method's top score decile who are also in first-ancestry PRS-CS's top decile; <i>r</i><sup>2</sup> is the squared correlation between scores and true genetic values.", common), t_nonnull),
  sec("Table 4. PRS-MIX against the seed rerun and the single panels", paste("Replicate-paired comparisons at the generating weight. Overlap and correlation are differences from the seed rerun (both first measured against first-ancestry PRS-CS); negative values mean PRS-MIX changes scores more than a new sampler seed alone. Null score variance is PRS-MIX relative to the seed rerun (geometric-mean percentage change). The last two columns are PRS-MIX's <i>r</i><sup>2</sup> minus single-panel PRS-CS's; positive values favour PRS-MIX. pp, percentage points.", common), t_contrast),
  "<section>", "<h2>Figure. Per-replicate differences across the sweep</h2>",
  sprintf("<p class=\"cap\">Each point is one replicate (simulated GWAS), with a violin for the distribution and a box for the median and quartiles; dashed line at zero. <b>a</b>, PRS-MIX minus the seed rerun in top-10%% overlap with first-ancestry PRS-CS (negative: PRS-MIX changes the top decile more than a new sampler seed). <b>b</b>, PRS-MIX minus first-ancestry PRS-CS in <i>r</i><sup>2</sup> with the true genetic value (positive: PRS-MIX more accurate). At most %d replicates per setting (see the tables for each count).</p>", max_rep),
  sprintf("<figure><img alt=\"Violin plots of per-replicate differences across the sweep\" src=\"data:image/png;base64,%s\"></figure>", fig_b64),
  "</section>",
  sprintf("<footer>Generated from <code>%s</code> by <code>inst/scripts/null_sim/13_sweep_summary.R</code> (at most %d replicates per setting; the Replicates column gives each setting's count) at PRS-MIX %s. Every number is in <code>sweep_tables.tsv</code>.</footer>",
          root, max_rep, substr(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), 1, 7)),
  "</main>", "</body>", "</html>")
writeLines(html, file.path(out_dir, "sweep_tables.html"))
message(sprintf("Wrote %s: %d settings, %d cells, replicates 1-%d", out_dir, length(settings), nrow(cells), max_rep))
