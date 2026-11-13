#!/usr/bin/env Rscript
# Journal-format tables for the 250 kb ancestry-pair x GWAS-composition sweep
# (4 pairs x 3 compositions), in the paper/<study>/ layout of paper/README.md.
#
# Rows are setting (pair, composition) x target ancestry. Panel A is the first
# ancestry of the pair (the one PRS-MIX's weight w refers to), panel B the
# second. PRS-MIX is run at the generating composition (w = share of A).
#
#   table_sweep_accuracy.tex  Table 1. r^2 with true genetic value and MSE for PRS-CS
#                             with panel A, PRS-MIX and PRS-CS with panel B,
#                             and PRS-MIX minus each single panel.
#   table_sweep_overlap.tex   Table 2. Top-10% overlap of PRS-MIX with each
#                             single panel, the seed-rerun baseline (the same
#                             single panel rerun with another sampler seed),
#                             and their difference.
#   table_sweep_cutoffs.tex   Table 3. That difference for the top 5%, 10%
#                             and 20% of target individuals.
#   table_sweep_settings.tex  Table 4. Simulation settings.
#   tables.html               The same tables for reading in a browser.
#   tables.tsv                Every number in the tables, unrounded, with n.
#   inputs/                   Frozen copies of the metric tables used.
#   manifest.json             Program, command, commits, input and output sha256.
#
# Cells are means over replicates with Monte Carlo standard errors (SD / sqrt(n));
# differences are paired within replicate (same simulated GWAS and target
# individuals), as in 12_paper_tables.R.
#
# Usage (from the PRS-MIX root):
#   Rscript inst/scripts/null_sim/18_sweep_paper_tables.R <sweep_dir> <out_dir> [reps=all] [phi_mode=auto] [provisional=no]
#
#   <sweep_dir>  e.g. inst/scripts/results/sim250_sweep_100rep, holding
#                <pair>/w<w>/metrics/nonnull_metrics.tsv for every setting.
#   reps         all, or N to use replicates 1..N of every setting.
#
# Writing under paper/ requires committed inputs and program, and at least 100
# replicates in every setting (the manuscript never quotes 20-replicate
# values). The one exception is a placeholder study: provisional=yes and a
# study folder whose name contains "provisional". Its tables carry a visible
# "Provisional placeholder" line and its manifest says provisional: true.
# Other output directories, such as a provisional 20-replicate preview under
# inst/scripts/results/, record the commit state in the manifest instead.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
args <- ns_args("Rscript 18_sweep_paper_tables.R <sweep_dir> <out_dir> [reps=all] [phi_mode=auto] [provisional=no]", 2)
sweep_dir <- args[[1]]
out_dir <- args[[2]]
reps_arg <- ns_arg(args, 3, "all")
phi_mode <- ns_arg(args, 4, "auto")
provisional <- ns_arg(args, 5, "no") == "yes"
program <- "inst/scripts/null_sim/18_sweep_paper_tables.R"
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
to_paper <- grepl("(^|/)paper/", normalizePath(out_dir))

sha <- function(path) sub("\\s.*$", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
git <- function(...) system2("git", c(...), stdout = TRUE)

pairs <- list(eur_asn = c("EUR", "ASN"), eur_afr = c("EUR", "AFR"), eur_sas = c("EUR", "SAS"), asn_afr = c("ASN", "AFR"))
ws <- c(0.9, 0.7, 0.5)
settings <- do.call(rbind, lapply(names(pairs), function(p) data.frame(pair = p, w = ws)))
settings$A <- vapply(settings$pair, function(p) pairs[[p]][1], "")
settings$B <- vapply(settings$pair, function(p) pairs[[p]][2], "")
settings$label <- sprintf("%s + %s, %.0f:%.0f", settings$A, settings$B, 100 * settings$w, 100 * (1 - settings$w))

# ---- Freeze inputs -------------------------------------------------------
inputs <- lapply(seq_len(nrow(settings)), function(i) {
  s <- settings[i, ]
  src <- file.path(sweep_dir, s$pair, sprintf("w%.1f", s$w), "metrics", "nonnull_metrics.tsv")
  dst <- file.path(out_dir, "inputs", sprintf("%s_w%.1f_nonnull_metrics.tsv", s$pair, s$w))
  file.copy(src, dst, overwrite = TRUE)
  stopifnot(sha(src) == sha(dst))
  list(path = dst, source = src, sha256 = sha(dst),
       source_commit = git("log", "-1", "--format=%H", "--", src),
       source_uncommitted = length(git("status", "--porcelain", "--", src)) > 0)
})
inputs_dirty <- any(vapply(inputs, `[[`, TRUE, "source_uncommitted"))
prog_dirty <- length(git("status", "--porcelain", "--", program, "inst/scripts/null_sim/common.R")) > 0
if (to_paper && (inputs_dirty || prog_dirty)) {
  stop("Writing under paper/ needs committed inputs and program; commit them first.", call. = FALSE)
}

# ---- Statistics ----------------------------------------------------------
cells <- list()
# scale = "log_pct": x holds per-replicate log ratios; report 100 (exp(mean) - 1)
# with SE 100 exp(mean) SE(mean log ratio), as in 12_paper_tables.R.
put <- function(table, i, target, column, x, note = "", scale = "mean") {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
  if (scale == "log_pct" && n) {
    se <- 100 * exp(m) * se
    m <- 100 * (exp(m) - 1)
  }
  cells[[length(cells) + 1]] <<- data.frame(
    table = table, setting = settings$label[i], target = target, column = column, n = n,
    mean = m, se = se, note = note, scale = scale)
}
n_reps <- integer(0)
for (i in seq_len(nrow(settings))) {
  s <- settings[i, ]
  nn <- ns_read_tsv(inputs[[i]]$path)
  nn <- nn[nn$phi_mode == phi_mode, ]
  if (reps_arg != "all") nn <- nn[nn$rep <= as.integer(reps_arg), ]
  mix <- sprintf("mix_w%.1f", s$w)
  sa <- paste0("single_", s$A); sb <- paste0("single_", s$B)
  arms <- c(mix, sa, sb, paste0(sa, "_seed2"), paste0(sb, "_seed2"))
  reps <- Reduce(intersect, lapply(arms, function(a) nn$rep[nn$arm == a]))
  n_reps[[s$label]] <- length(reps)
  v <- function(arm, col) { y <- nn[nn$arm == arm, ]; y[[col]][match(reps, y$rep)] }
  for (tg in c(s$A, s$B)) {
    r2 <- function(a) v(a, paste0("r2_", tg))
    put("accuracy", i, tg, "r2_A", r2(sa))
    put("accuracy", i, tg, "r2_mix", r2(mix))
    put("accuracy", i, tg, "r2_B", r2(sb))
    put("accuracy", i, tg, "d_r2_A", r2(mix) - r2(sa))
    put("accuracy", i, tg, "d_r2_B", r2(mix) - r2(sb))
    lmse <- function(a) log(v(mix, paste0("mse_gv_", tg)) / v(a, paste0("mse_gv_", tg)))
    put("accuracy", i, tg, "mse_pct_A", lmse(sa), scale = "log_pct")
    put("accuracy", i, tg, "mse_pct_B", lmse(sb), scale = "log_pct")
    for (ref in c("A", "B")) {
      rp <- s[[ref]]
      for (cut in c("05", "10", "20")) {
        col <- sprintf("top%s_%s_%s", cut, rp, tg)
        ov <- 100 * v(mix, col)
        base <- 100 * v(paste0("single_", rp, "_seed2"), col)
        if (cut == "10") {
          put("overlap", i, tg, paste0("ov_", ref), ov)
          put("overlap", i, tg, paste0("base_", ref), base)
          put("overlap", i, tg, paste0("gap_", ref), ov - base)
        }
        put("cutoffs", i, tg, sprintf("gap%s_%s", cut, ref), ov - base)
      }
    }
  }
}
cells <- do.call(rbind, cells)
placeholder <- provisional && grepl("provisional", basename(normalizePath(out_dir)))
if (provisional && !placeholder) stop("provisional=yes needs a study folder whose name contains 'provisional'.", call. = FALSE)
if (to_paper && any(n_reps < 100) && !placeholder) {
  stop(sprintf("Writing under paper/ needs >= 100 replicates per setting; fewest is %d (%s).",
               min(n_reps), names(which.min(n_reps))), call. = FALSE)
}
ns_write_tsv(cells, file.path(out_dir, "tables.tsv"))
rep_text <- if (length(unique(n_reps)) == 1) as.character(n_reps[[1]]) else sprintf("%d to %d", min(n_reps), max(n_reps))

# ---- Table definitions (shared by LaTeX and HTML) ------------------------
# Rows are grouped by ancestry pair. Each pair's block opens with a header row
# naming its real panels (e.g. "EUR panel"), so no A/B code is needed. Text
# helpers take tex = TRUE/FALSE. pa/pb are the pair's ancestries.
col <- function(id, digits, signed = FALSE) list(id = id, digits = digits, signed = signed)
tx <- function(tex, t, h) if (tex) t else h
minus <- function(tex) tx(tex, "$-$", "&minus;")
tables <- list(
  accuracy = list(
    cols = list(col("r2_A", 3), col("r2_mix", 3), col("r2_B", 3), col("d_r2_A", 3, TRUE), col("d_r2_B", 3, TRUE),
                col("mse_pct_A", 1, TRUE), col("mse_pct_B", 1, TRUE)),
    groups = function(tex) list(
      list(tx(tex, "$r^2$ with true genetic value (higher is better)", "<i>r</i><sup>2</sup> with true genetic value (higher is better)"), 3),
      list(tx(tex, "$r^2$, PRS-MIX minus single panel (positive favours PRS-MIX)", "<i>r</i><sup>2</sup>, PRS-MIX minus single panel (positive favours PRS-MIX)"), 2),
      list(tx(tex, "MSE, PRS-MIX vs.\\ single panel (\\%; negative favours PRS-MIX)", "MSE, PRS-MIX vs. single panel (%; negative favours PRS-MIX)"), 2)),
    heads = NULL,
    block = function(pa, pb, tex) list(list(paste(pa, "panel"), 1), list("PRS-MIX", 1), list(paste(pb, "panel"), 1),
                                       list(tx(tex, paste0("vs.\\ ", pa, " panel"), paste("vs.", pa, "panel")), 1),
                                       list(tx(tex, paste0("vs.\\ ", pb, " panel"), paste("vs.", pb, "panel")), 1),
                                       list(tx(tex, paste0("vs.\\ ", pa, " panel"), paste("vs.", pa, "panel")), 1),
                                       list(tx(tex, paste0("vs.\\ ", pb, " panel"), paste("vs.", pb, "panel")), 1))),
  overlap = list(
    cols = list(col("ov_A", 1), col("base_A", 1), col("gap_A", 1, TRUE), col("ov_B", 1), col("base_B", 1), col("gap_B", 1, TRUE)),
    groups = function(tex) rep(list(list(tx(tex, "Top-10\\% overlap with the panel (\\%)", "Top-10% overlap with the panel (%)"), 2),
                                    list("Difference", 1)), 2),
    heads = function(tex) rep(c("PRS-MIX", "Rerun, new seed", "Percentage points"), 2),
    block = function(pa, pb, tex) list(list(paste("Compared with the", pa, "panel"), 3), list(paste("Compared with the", pb, "panel"), 3))),
  cutoffs = list(
    cols = list(col("gap05_A", 1, TRUE), col("gap10_A", 1, TRUE), col("gap20_A", 1, TRUE),
                col("gap05_B", 1, TRUE), col("gap10_B", 1, TRUE), col("gap20_B", 1, TRUE)),
    groups = function(tex) rep(list(list("Overlap difference, PRS-MIX minus rerun (percentage points)", 3)), 2),
    heads = function(tex) rep(tx(tex, c("Top 5\\%", "Top 10\\%", "Top 20\\%"), c("Top 5%", "Top 10%", "Top 20%")), 2),
    block = function(pa, pb, tex) list(list(paste("Compared with the", pa, "panel"), 3), list(paste("Compared with the", pb, "panel"), 3)))
)
fmt_num <- function(x, d, signed, tex) {
  s <- formatC(abs(x), format = "f", digits = d)
  if (as.numeric(s) == 0) return(s)  # no sign on a value that rounds to zero
  sign <- if (x < 0) minus(tex) else if (signed && x > 0) "+" else ""
  paste0(sign, s)
}
cell_text <- function(table, label, target, c, tex) {
  r <- cells[cells$table == table & cells$setting == label & cells$target == target & cells$column == c$id, ]
  if (nrow(r) == 0 || is.na(r$mean)) return(tx(tex, "--", "&ndash;"))
  if (is.na(r$se)) return(fmt_num(r$mean, c$digits, c$signed, tex))
  unit <- 10^-c$digits
  se <- if (r$se < unit / 2) paste0(tx(tex, "$<$", "&lt;"), formatC(unit, format = "f", digits = c$digits))
        else formatC(r$se, format = "f", digits = c$digits)
  sprintf("%s (%s)", fmt_num(r$mean, c$digits, c$signed, tex), se)
}
mix_label <- function(i, tex) {
  pct <- tx(tex, "\\%", "%")
  sprintf("%.0f%s %s, %.0f%s %s", 100 * settings$w[i], pct, settings$A[i], 100 * (1 - settings$w[i]), pct, settings$B[i])
}
pair_label <- function(p) paste(pairs[[p]][1], "and", pairs[[p]][2])

# ---- LaTeX ---------------------------------------------------------------
latex_table <- function(name) {
  t <- tables[[name]]
  nc <- length(t$cols)
  spans <- function(items) {  # multicolumn cells and cmidrules for a list of (text, span)
    l <- vapply(items, `[[`, 1, 2)
    ends <- cumsum(l) + 2
    starts <- ends - l + 1
    list(cells = vapply(items, function(x) if (x[[2]] == 1) x[[1]] else sprintf("\\multicolumn{%d}{c}{%s}", x[[2]], x[[1]]), ""),
         rules = paste(sprintf("\\cmidrule(lr){%d-%d}", starts, ends), collapse = ""))
  }
  g <- spans(t$groups(TRUE))
  head <- c(paste(paste(c("GWAS ancestry mix", "Scored in", g$cells), collapse = " & "), "\\\\"), g$rules,
            if (!is.null(t$heads)) paste(paste(c("", "", t$heads(TRUE)), collapse = " & "), "\\\\"))
  body <- unlist(lapply(names(pairs), function(p) {
    pa <- pairs[[p]][1]; pb <- pairs[[p]][2]
    b <- spans(t$block(pa, pb, TRUE))
    rows <- unlist(lapply(which(settings$pair == p), function(i) {
      c(sprintf("%s & %s & %s \\\\", mix_label(i, TRUE), pa,
                paste(vapply(t$cols, function(c) cell_text(name, settings$label[i], pa, c, TRUE), ""), collapse = " & ")),
        sprintf(" & %s & %s \\\\", pb,
                paste(vapply(t$cols, function(c) cell_text(name, settings$label[i], pb, c, TRUE), ""), collapse = " & ")),
        if (i != max(which(settings$pair == p))) "\\addlinespace[2pt]")
    }))
    c(if (p != names(pairs)[1]) "\\addlinespace", sprintf("\\multicolumn{2}{l}{\\textbf{%s}} & %s \\\\", pair_label(p), paste(b$cells, collapse = " & ")),
      b$rules, rows)
  }))
  c("% Generated by inst/scripts/null_sim/18_sweep_paper_tables.R; do not edit.",
    sprintf("\\begin{tabular}{ll%s}", strrep("r", nc)), "\\toprule",
    head, "\\midrule", body, "\\bottomrule",
    if (any(n_reps < 100)) sprintf("\\multicolumn{%d}{l}{\\footnotesize\\textit{Provisional placeholder: %s replicates per setting; not for quoting.}} \\\\", nc + 2, rep_text),
    "\\end{tabular}")
}
settings_rows <- list(
  c("Ancestry pairs", "EUR and ASN, EUR and AFR, EUR and SAS, ASN and AFR"),
  c("GWAS ancestry mixes", "90:10, 70:30 and 50:50 for each pair"),
  c("Chromosome", "22, HapMap3 SNPs"),
  c("LD blocks", "fixed 250 kb, shared by both panels"),
  c("Simulated GWAS", "1\\% causal SNPs with effects shared across ancestries; nominal $N$ = 200,000"),
  c("Replicates", sprintf("%s per setting", rep_text)),
  c("PRS-MIX blend", "each panel weighted by its ancestry's share of the GWAS"),
  c("PRS-CS", sprintf("$\\phi$ %s, sampler defaults (1,000 iterations, 500 burn-in, thinning 5)", phi_mode)),
  c("Sampler seeds", "1000 + replicate; rerun 2000 + replicate"))
latex_settings <- c("% Generated by inst/scripts/null_sim/18_sweep_paper_tables.R; do not edit.",
                    "\\begin{tabular}{ll}", "\\toprule", "Item & Value \\\\", "\\midrule",
                    vapply(settings_rows, function(r) sprintf("%s & %s \\\\", r[1], r[2]), ""),
                    "\\bottomrule", "\\end{tabular}")
writeLines(latex_table("accuracy"), file.path(out_dir, "table_sweep_accuracy.tex"))
writeLines(latex_table("overlap"), file.path(out_dir, "table_sweep_overlap.tex"))
writeLines(latex_table("cutoffs"), file.path(out_dir, "table_sweep_cutoffs.tex"))
writeLines(latex_settings, file.path(out_dir, "table_sweep_settings.tex"))

# ---- HTML ----------------------------------------------------------------
tex2html <- function(s) {
  s <- gsub("\\\\%", "%", s)
  s <- gsub("\\$\\\\phi\\$", "<i>&phi;</i>", s)
  s <- gsub("\\$N\\$", "<i>N</i>", s)
  s
}
html_table <- function(name) {
  t <- tables[[name]]
  nc <- length(t$cols)
  th_span <- function(items, cls) paste(vapply(items, function(x)
    sprintf("<th%s%s>%s</th>", if (x[[2]] > 1) sprintf(" colspan=\"%d\"", x[[2]]) else "", cls, x[[1]]), ""), collapse = "")
  rs <- if (is.null(t$heads)) "" else " rowspan=\"2\""
  head <- c(sprintf("<tr class=\"top\"><th%s>GWAS ancestry mix</th><th%s>Scored in</th>%s</tr>", rs, rs, th_span(t$groups(FALSE), " class=\"grp\"")),
            if (!is.null(t$heads)) sprintf("<tr class=\"sub\">%s</tr>", paste(sprintf("<th>%s</th>", t$heads(FALSE)), collapse = "")))
  body <- unlist(lapply(names(pairs), function(p) {
    pa <- pairs[[p]][1]; pb <- pairs[[p]][2]
    rows <- unlist(lapply(which(settings$pair == p), function(i) {
      c(sprintf("<tr class=\"mixrow%s\"><td>%s</td><td>%s</td>%s</tr>", if (i == min(which(settings$pair == p))) " first" else "", mix_label(i, FALSE), pa,
                paste(sprintf("<td>%s</td>", vapply(t$cols, function(c) cell_text(name, settings$label[i], pa, c, FALSE), "")), collapse = "")),
        sprintf("<tr><td></td><td>%s</td>%s</tr>", pb,
                paste(sprintf("<td>%s</td>", vapply(t$cols, function(c) cell_text(name, settings$label[i], pb, c, FALSE), "")), collapse = "")))
    }))
    c("<tbody>", sprintf("<tr class=\"block\"><th colspan=\"2\" scope=\"rowgroup\">%s</th>%s</tr>", pair_label(p), th_span(t$block(pa, pb, FALSE), "")),
      rows, "</tbody>")
  }))
  c("<table>", "<thead>", head, "</thead>", body, "</table>")
}
html_settings <- c("<table class=\"kv\">", "<thead>", "<tr class=\"top\"><th>Item</th><th>Value</th></tr>", "</thead>", "<tbody>",
                   vapply(settings_rows, function(r) sprintf("<tr><td>%s</td><td>%s</td></tr>", r[1], tex2html(r[2])), ""),
                   "</tbody>", "</table>")
glossary <- c(
  "<section class=\"howto\">", "<h2>How to read these tables</h2>", "<dl>",
  "<dt>EUR panel, ASN panel, &hellip;</dt><dd>Standard PRS-CS run with one ancestry's LD reference panel.</dd>",
  "<dt>PRS-MIX</dt><dd>PRS-CS run with an LD panel blended from the pair's two panels, each weighted by its ancestry's share of the GWAS (for example 70% EUR panel and 30% ASN panel for a 70% EUR, 30% ASN GWAS).</dd>",
  "<dt>GWAS ancestry mix</dt><dd>The ancestry make-up of the simulated GWAS. Each pair is tested at three mixes.</dd>",
  "<dt>Scored in</dt><dd>The held-out individuals whose scores are evaluated. They are in neither the GWAS nor the LD panels.</dd>",
  "<dt><i>r</i><sup>2</sup> with true genetic value</dt><dd>Squared correlation between people's scores and their true simulated genetic values. 1 is perfect, 0 is no information.</dd>",
  "<dt>MSE</dt><dd>Mean squared error between people's scores and their true genetic values. Unlike <i>r</i><sup>2</sup>, it also penalises scores that are on the wrong scale, for example spread far more widely than the true values. The table gives PRS-MIX's MSE as a percentage change from each single panel's: &minus;40% means PRS-MIX's error is 40% smaller.</dd>",
  "<dt>Top-10% overlap</dt><dd>Of the people PRS-MIX ranks in its top 10%, the share that the single panel also ranks in its top 10%.</dd>",
  "<dt>Rerun, new seed</dt><dd>The same single panel run again with a different random seed. Its overlap with the original run shows how much the top 10% changes from sampler randomness alone. It is the yardstick for PRS-MIX.</dd>",
  "<dt>Percentage points</dt><dd>The plain difference between two percentages: 88% minus 91% is &minus;3 percentage points. Negative means PRS-MIX changes who is in the top 10% more than rerunning the same panel does.</dd>",
  "<dt>Numbers in parentheses</dt><dd>Standard errors over the simulation replicates. Differences are paired: within a replicate, every method analyses the same GWAS and scores the same people.</dd>",
  "</dl>", "</section>")
common_caption <- sprintf(paste(
  "Chromosome 22, fixed 250 kb LD blocks. Means over %s replicates, standard errors in parentheses; differences are",
  "paired within replicate."), rep_text)
captions <- list(
  accuracy = paste("How well each method's scores track true genetic value. <i>r</i><sup>2</sup> measures how well scores rank",
                   "people; MSE also checks that scores are on the right scale. The <i>r</i><sup>2</sup> differences subtract each",
                   "single panel's <i>r</i><sup>2</sup> from PRS-MIX's. The MSE columns give PRS-MIX's MSE as a percentage change from",
                   "the single panel's, averaged on the log scale over replicates (geometric mean) with delta-method standard errors.",
                   common_caption),
  overlap = paste("For each single panel: how much of PRS-MIX's top 10% that panel also puts in its top 10%, next to the same",
                  "overlap for a rerun of the panel with a new seed. The difference shows how far PRS-MIX moves the top 10%",
                  "beyond sampler randomness; more negative means a bigger change.", common_caption),
  cutoffs = paste("The difference from Table 2 (PRS-MIX overlap minus rerun overlap, in percentage points) for the top 5%,",
                  "10% and 20% of scores.", common_caption),
  settings = "Simulation and analysis settings for the sweep.")
titles <- c(accuracy = "Table 1. Accuracy: how well scores track true genetic value",
            overlap = "Table 2. Who is in the top 10%: PRS-MIX vs. each single panel, against rerun randomness",
            cutoffs = "Table 3. The same comparison at the top 5%, 10% and 20%",
            settings = "Table 4. Simulation settings")
section <- function(name, tab) c("<section>", sprintf("<h2>%s</h2>", titles[[name]]),
                                 sprintf("<p class=\"cap\">%s</p>", captions[[name]]), "<div class=\"wrap\">", tab, "</div>", "</section>")
head_commit <- git("rev-parse", "HEAD")
status <- if (inputs_dirty || prog_dirty) sprintf(" with uncommitted %s", paste(c(if (prog_dirty) "program", if (inputs_dirty) "inputs"), collapse = " and ")) else ""
prov_html <- if (any(n_reps < 100)) sprintf("<p class=\"prov\">Provisional: %s replicates per setting. Manuscript numbers come from the 100-replicate runs.</p>", rep_text) else ""
sans <- "-apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif"
html <- c(
  "<!doctype html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
  sprintf("<title>%s tables</title>", basename(normalizePath(out_dir))),
  "<!-- Generated by inst/scripts/null_sim/18_sweep_paper_tables.R; do not edit. -->",
  "<style>",
  "body { margin: 0; padding: 40px 24px 56px; background: #fff; color: #111; font: 15px/1.45 Georgia, \"Times New Roman\", Times, serif; }",
  "main { max-width: 1180px; margin: 0 auto; display: flex; flex-direction: column; gap: 44px; }",
  sprintf("h2 { font: 600 17px/1.3 %s; margin: 0 0 8px; }", sans),
  ".cap { margin: 0 0 14px; font-size: 14px; max-width: none; }",
  sprintf(".prov { margin: 0; font: 600 13px/1.4 %s; color: #8a3020; }", sans),
  ".howto { border-top: 2px solid #111; border-bottom: 1px solid #111; padding: 14px 0 6px; }",
  ".howto dl { display: grid; grid-template-columns: max-content 1fr; gap: 6px 20px; margin: 0; font-size: 14px; }",
  ".howto dt { font-weight: 700; } .howto dd { margin: 0; }",
  "@media (max-width: 640px) { .howto dl { grid-template-columns: 1fr; } .howto dd { margin-bottom: 6px; } }",
  ".wrap { overflow-x: auto; }",
  "table { border-collapse: collapse; width: 100%; font-variant-numeric: lining-nums tabular-nums; font-size: 14.5px; }",
  "th, td { padding: 4px 10px; text-align: center; white-space: nowrap; vertical-align: baseline; }",
  "th:nth-child(-n+2), td:nth-child(-n+2) { text-align: left; }",
  "table.kv td:last-child { white-space: normal; }",
  "thead th { font-weight: 700; vertical-align: bottom; }",
  "thead tr.top th { border-top: 2px solid #111; padding-top: 6px; }",
  "thead tr:last-child th, thead tr.top th[rowspan] { border-bottom: 1px solid #111; padding-bottom: 6px; }",
  "thead th.grp { border-bottom: 1px solid #111; }",
  "tr.block th { font-weight: 700; padding-top: 14px; border-bottom: 1px solid #999; text-align: center; }",
  "tr.block th:first-child { text-align: left; }",
  "tbody:first-of-type tr.block th { padding-top: 8px; }",
  "tr.mixrow td { padding-top: 9px; }",
  "tr.mixrow.first td { padding-top: 6px; }",
  "table:not(.kv) tbody:last-of-type tr:last-child td, table.kv tbody tr:last-child td { border-bottom: 2px solid #111; padding-bottom: 6px; }",
  sprintf("footer { font: 12px/1.5 %s; color: #666; }", sans),
  "code { font-size: 11.5px; }",
  "</style>", "</head>", "<body>", "<main>", prov_html, glossary,
  section("accuracy", html_table("accuracy")),
  section("overlap", html_table("overlap")),
  section("cutoffs", html_table("cutoffs")),
  section("settings", html_settings),
  sprintf(paste("<footer>Rendered from <code>%s/inputs/</code> by <code>%s</code> at PRS-MIX %s%s. The <code>.tex</code> files in",
                "the same folder hold the same values, and <code>tables.tsv</code> holds them unrounded with replicate",
                "counts.</footer>"), out_dir, program, substr(head_commit, 1, 7), status),
  "</main>", "</body>", "</html>")
writeLines(html, file.path(out_dir, "tables.html"))

# ---- Manifest ------------------------------------------------------------
out_files <- c(table_sweep_accuracy = "table_sweep_accuracy.tex", table_sweep_overlap = "table_sweep_overlap.tex",
               table_sweep_cutoffs = "table_sweep_cutoffs.tex", table_sweep_settings = "table_sweep_settings.tex",
               tables_tsv = "tables.tsv", tables_html = "tables.html")
outputs <- lapply(out_files, function(f) list(path = file.path(out_dir, f), sha256 = sha(file.path(out_dir, f))))
manifest <- list(
  study = basename(normalizePath(out_dir)), program = program,
  command = paste("Rscript", program, paste(shQuote(args), collapse = " ")),
  asset_commit = head_commit, program_uncommitted = prog_dirty, inputs_uncommitted = inputs_dirty,
  provisional = any(n_reps < 100),
  experiment_commit = unique(vapply(inputs, `[[`, "", "source_commit")),
  r_version = R.version.string,
  inputs = unname(lapply(inputs, function(x) x[c("path", "source", "source_commit", "sha256")])),
  outputs = outputs,
  settings = list(
    chromosome = 22, block_layout = "fixed 250 kb", pairs = unname(lapply(pairs, as.list)), compositions_A = ws,
    replicates_used = as.list(n_reps), replicate_filter = reps_arg, phi_mode = phi_mode,
    prs_mix_weight = "generating composition (share of panel A's ancestry)",
    seed_control = "single-panel arm rerun with sampler seed 2000 + rep (main runs use 1000 + rep)",
    uncertainty = "mean over replicates with Monte Carlo SE (SD / sqrt(n)); paired differences use the SE of the per-replicate difference",
    metric_definitions = list(
      r2 = "squared correlation between target scores and true genetic values, non-null GWAS",
      mse_pct = "100 * (exp(mean over replicates of log(mse_gv of PRS-MIX / mse_gv of the single panel)) - 1); mse_gv = mean of ((s - mean s) - (g - mean g))^2 (07_metrics_nonnull.R); delta-method SE",
      overlap = "share of target individuals in PRS-MIX's top score fraction also in the single panel's top fraction (top05/top10/top20 columns of 07_metrics_nonnull.R)",
      seed_baseline = "the same overlap between the single panel's seed rerun and the single panel",
      difference = "per-replicate PRS-MIX overlap minus seed-rerun overlap, percentage points")),
  notes = paste("Replicates missing any of the five arms used are dropped from that setting; n per cell is in tables.tsv.",
                "Experiment runs were made on RedHawk with the null-simulation job scripts (cluster/jobs/null_sim, not version-controlled)."))
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(sprintf("Wrote %s (%d cells, %s replicates) at PRS-MIX %s%s", out_dir, nrow(cells), rep_text, substr(head_commit, 1, 7), status))
