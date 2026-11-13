#!/usr/bin/env Rscript
# Curated manuscript tables for one null/non-null simulation study (plan
# section 7), written into paper/<study>/ for import by the manuscript
# repository.
#
# Methods (rows): EUR-only PRS-CS (reference), the same run with another sampler
# seed (run-to-run control), PRS-MIX at w = 0.9, ..., 0.1 (weight on pops[1]),
# and ASN-only PRS-CS. PRS-MIX at w = 1.0 and 0.0 is identical to the
# single-panel runs, so those weights are not repeated.
#
#   table_null.tex           Table 1. Null GWAS: excess LD-fit loss over the
#                            generating weight; null score variance relative to
#                            EUR-only PRS-CS, by target ancestry.
#   table_nonnull.tex        Table 2. Non-null GWAS: top-decile overlap and score
#                            correlation with EUR-only PRS-CS; r^2 with the true
#                            genetic value; by target ancestry.
#   table_seed_contrast.tex  Table 3. Replicate-paired differences, each method
#                            minus the EUR-only seed rerun.
#   table_settings.tex       Table 4. Simulation settings.
#   tables.html              The same four tables for reading in a browser.
#   tables.tsv               Every number in the tables, unrounded, with n.
#   inputs/                  Frozen copies of the metric tables used.
#   manifest.json            Program, command, commits, input and output sha256.
#
# Cells are means over replicates with Monte Carlo standard errors (SD / sqrt(n)).
# Comparisons are paired within replicate (same simulated GWAS and target
# individuals); their SE is the SE of the paired difference. Score-variance
# ratios are averaged on the log scale and reported as a percentage change,
# 100 (exp(mean log ratio) - 1), with a delta-method SE; averaging the ratio
# itself would bias it upward (two equally variable runs would average > 0).
# The reference row reads 0 (ref.), 100 (ref.) or 1 (ref.): its value by definition.
#
# Usage (from the PRS-MIX root, with the inputs committed):
#   Rscript inst/scripts/null_sim/12_paper_tables.R <metrics_dir> <out_dir> \
#     [phi_mode=auto] [w0=0.7] [pops=EUR,ASN]

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
args <- ns_args("Rscript 12_paper_tables.R <metrics_dir> <out_dir> [phi_mode=auto] [w0=0.7] [pops=EUR,ASN]", 2)
metrics_dir <- args[[1]]
out_dir <- args[[2]]
phi_mode <- ns_arg(args, 3, "auto")
w0 <- as.numeric(ns_arg(args, 4, "0.7"))
pops <- toupper(strsplit(ns_arg(args, 5, "EUR,ASN"), ",")[[1]])
stopifnot(length(pops) == 2)
pa <- pops[1]
pb <- pops[2]
program <- "inst/scripts/null_sim/12_paper_tables.R"

sha <- function(path) sub("\\s.*$", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
git <- function(...) system2("git", c(...), stdout = TRUE)

# ---- Freeze inputs -------------------------------------------------------
input_names <- c("prscs_metrics.tsv", "nonnull_metrics.tsv", "ld_fit.tsv")
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
inputs <- lapply(input_names, function(f) {
  src <- file.path(metrics_dir, f)
  dst <- file.path(out_dir, "inputs", f)
  file.copy(src, dst, overwrite = TRUE)
  stopifnot(sha(src) == sha(dst))
  list(path = dst, source = src, sha256 = sha(dst),
       source_commit = git("log", "-1", "--format=%H", "--", src),
       source_uncommitted = length(git("status", "--porcelain", "--", src)) > 0)
})
names(inputs) <- input_names
if (any(vapply(inputs, `[[`, TRUE, "source_uncommitted"))) {
  stop("Input metric tables have uncommitted changes; commit them first.", call. = FALSE)
}
pm <- ns_read_tsv(inputs[["prscs_metrics.tsv"]]$path)
nn <- ns_read_tsv(inputs[["nonnull_metrics.tsv"]]$path)
ld <- ns_read_tsv(inputs[["ld_fit.tsv"]]$path)
pm <- pm[pm$phi_mode == phi_mode, ]
nn <- nn[nn$phi_mode == phi_mode, ]
ld <- ld[abs(ld$delta - 0.001) < 1e-12, ]

# ---- Statistics ----------------------------------------------------------
ref <- paste0("single_", pa)
seed_arm <- paste0(ref, "_seed2")
mix_arms <- sprintf("mix_w%.1f", seq(0.9, 0.1, by = -0.1))
methods <- c(ref, seed_arm, mix_arms, paste0("single_", pb))
weight_of <- function(arm) ns_arm_weight(arm, pops)
setting <- sprintf("%s + %s, %.0f:%.0f", pa, pb, 100 * w0, 100 * (1 - w0))

# Per-replicate values of `col` for `arm`, named by replicate.
per_rep <- function(tab, arm, col) {
  x <- tab[tab$arm == arm, ]
  stats::setNames(x[[col]], x$rep)
}
paired <- function(a, b) {
  r <- intersect(names(a), names(b))
  a[r] - b[r]
}
ld_excess <- function(arm) {
  a <- ld[abs(ld$w - weight_of(arm)) < 1e-9, ]
  b <- ld[abs(ld$w - w0) < 1e-9, ]
  paired(stats::setNames(a$L, a$rep), stats::setNames(b$L, b$rep))
}
log_var_ratio <- function(arm, base, p) {  # per-replicate log(var_arm / var_base)
  col <- paste0("score_var_", p)
  paired(log(per_rep(pm, arm, col)), log(per_rep(pm, base, col)))
}

cells <- list()
# scale = "log_pct": x holds per-replicate log ratios; report 100 (exp(mean) - 1)
# with SE 100 exp(mean) SE(mean log ratio).
put <- function(table, arm, column, x, note = "", scale = "mean") {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
  if (scale == "log_pct" && n) {
    se <- 100 * exp(m) * se
    m <- 100 * (exp(m) - 1)
  }
  cells[[length(cells) + 1]] <<- data.frame(
    table = table, setting = setting, arm = arm, column = column, n = n,
    mean = m, se = se, note = note, scale = scale)
}
for (arm in methods) {
  is_ref <- arm == ref
  is_seed <- arm == seed_arm
  # Table 1: null GWAS.
  put("null", arm, "replicates", numeric(0), as.character(sum(pm$arm == arm)))
  if (is_seed) put("null", arm, "ld_fit_excess", numeric(0), "n/a") else put("null", arm, "ld_fit_excess", ld_excess(arm))
  for (p in pops) {
    if (is_ref) put("null", arm, paste0("null_var_pct_", p), numeric(0), "0 (ref.)")
    else put("null", arm, paste0("null_var_pct_", p), log_var_ratio(arm, ref, p), scale = "log_pct")
  }
  # Table 2: non-null GWAS.
  put("nonnull", arm, "replicates", numeric(0), as.character(sum(nn$arm == arm)))
  for (p in pops) {
    if (is_ref) {
      put("nonnull", arm, paste0("top10_pct_", p), numeric(0), "100 (ref.)")
      put("nonnull", arm, paste0("cor_", p), numeric(0), "1 (ref.)")
    } else {
      put("nonnull", arm, paste0("top10_pct_", p), 100 * per_rep(nn, arm, sprintf("top10_%s_%s", pa, p)))
      put("nonnull", arm, paste0("cor_", p), per_rep(nn, arm, sprintf("cor_%s_%s", pa, p)))
    }
  }
  for (p in pops) put("nonnull", arm, paste0("r2_", p), per_rep(nn, arm, paste0("r2_", p)))
  # Table 3: method minus the seed rerun, paired by replicate.
  if (!is_ref && !is_seed) {
    d_var <- lapply(stats::setNames(pops, pops), function(p) log_var_ratio(arm, seed_arm, p))
    d_top <- lapply(stats::setNames(pops, pops), function(p) {
      col <- sprintf("top10_%s_%s", pa, p)
      100 * paired(per_rep(nn, arm, col), per_rep(nn, seed_arm, col))
    })
    d_cor <- lapply(stats::setNames(pops, pops), function(p) {
      col <- sprintf("cor_%s_%s", pa, p)
      1000 * paired(per_rep(nn, arm, col), per_rep(nn, seed_arm, col))
    })
    put("seed_contrast", arm, "pairs", numeric(0),
        sprintf("%d / %d", length(d_var[[1]]), length(d_top[[1]])))
    for (p in pops) put("seed_contrast", arm, paste0("d_null_var_", p), d_var[[p]], scale = "log_pct")
    for (p in pops) put("seed_contrast", arm, paste0("d_top10_", p), d_top[[p]])
    for (p in pops) put("seed_contrast", arm, paste0("d_cor_x1000_", p), d_cor[[p]])
  }
}
cells <- do.call(rbind, cells)
ns_write_tsv(cells, file.path(out_dir, "tables.tsv"))

# ---- Table definitions (shared by LaTeX and HTML) ------------------------
# Each column: id, digits, signed (show + for positive), header group, header.
col <- function(id, digits, group = "", head, signed = FALSE) list(id = id, digits = digits, group = group, head = head, signed = signed)
tables <- list(
  null = list(
    rows = methods, label_head = "Method",
    cols = list(col("replicates", 0, "", "Replicates"),
                col("ld_fit_excess", 3, "", sprintf("LD-fit loss above $w = %.1f$", w0)),
                col(paste0("null_var_pct_", pa), 1, "nullvar", sprintf("%s targets", pa), signed = TRUE),
                col(paste0("null_var_pct_", pb), 1, "nullvar", sprintf("%s targets", pb), signed = TRUE)),
    groups = list(nullvar = sprintf("Null score variance vs.\\ %s-only (\\%%)", pa))),
  nonnull = list(
    rows = methods, label_head = "Method",
    cols = list(col("replicates", 0, "", "Replicates"),
                col(paste0("top10_pct_", pa), 1, "top", sprintf("%s", pa)),
                col(paste0("top10_pct_", pb), 1, "top", sprintf("%s", pb)),
                col(paste0("cor_", pa), 3, "cor", sprintf("%s", pa)),
                col(paste0("cor_", pb), 3, "cor", sprintf("%s", pb)),
                col(paste0("r2_", pa), 3, "r2", sprintf("%s", pa)),
                col(paste0("r2_", pb), 3, "r2", sprintf("%s", pb))),
    groups = list(top = sprintf("Top-10\\%% overlap with %s-only (\\%%)", pa),
                  cor = sprintf("Correlation with %s-only", pa),
                  r2 = "$r^2$ with true genetic value")),
  seed_contrast = list(
    rows = c(mix_arms, paste0("single_", pb)), label_head = "Comparison",
    cols = list(col("pairs", 0, "", "Pairs (null / non-null)"),
                col(paste0("d_null_var_", pa), 1, "dvar", pa, signed = TRUE),
                col(paste0("d_null_var_", pb), 1, "dvar", pb, signed = TRUE),
                col(paste0("d_top10_", pa), 2, "dtop", pa, signed = TRUE),
                col(paste0("d_top10_", pb), 2, "dtop", pb, signed = TRUE),
                col(paste0("d_cor_x1000_", pa), 2, "dcor", pa, signed = TRUE),
                col(paste0("d_cor_x1000_", pb), 2, "dcor", pb, signed = TRUE)),
    groups = list(dvar = "Null score variance (\\%)", dtop = "Top-10\\% overlap (pp)",
                  dcor = "Correlation $\\times 10^3$"))
)
method_label <- function(arm, tex) {
  star <- function(w) if (abs(w - w0) < 1e-9) (if (tex) "$^{*}$" else "<sup>*</sup>") else ""
  wl <- function(w) if (tex) sprintf("PRS-MIX, $w = %.1f$%s", w, star(w)) else sprintf("PRS-MIX, <i>w</i> = %.1f%s", w, star(w))
  if (arm == ref) return(sprintf("PRS-CS, %s panel", pa))
  if (arm == seed_arm) return(sprintf("PRS-CS, %s panel, new seed", pa))
  if (grepl("^single_", arm)) return(sprintf("PRS-CS, %s panel", sub("^single_", "", arm)))
  wl(weight_of(arm))
}
contrast_label <- function(arm, tex) {
  minus <- if (tex) " $-$ " else " &minus; "
  paste0(method_label(arm, tex), minus, sprintf("%s rerun", pa))
}
fmt_num <- function(x, d, signed, tex) {
  s <- formatC(abs(x), format = "f", digits = d)
  if (as.numeric(s) == 0) return(s)  # no sign on a value that rounds to zero
  sign <- if (x < 0) (if (tex) "$-$" else "&minus;") else if (signed && x > 0) "+" else ""
  paste0(sign, s)
}
cell_text <- function(table, arm, c, tex) {
  r <- cells[cells$table == table & cells$arm == arm & cells$column == c$id, ]
  if (nrow(r) == 0) return(if (tex) "--" else "&ndash;")
  if (r$note == "n/a") return(if (tex) "--" else "&ndash;")
  if (nzchar(r$note)) return(r$note)
  if (is.na(r$se)) return(fmt_num(r$mean, c$digits, c$signed, tex))
  if (r$mean == 0 && r$se == 0) return("0")  # zero by definition (generating weight)
  unit <- 10^-c$digits
  se <- if (r$se < unit / 2) paste0(if (tex) "$<$" else "&lt;", formatC(unit, format = "f", digits = c$digits))
        else formatC(r$se, format = "f", digits = c$digits)
  sprintf("%s (%s)", fmt_num(r$mean, c$digits, c$signed, tex), se)
}
row_label <- function(name, arm, tex) if (name == "seed_contrast") contrast_label(arm, tex) else method_label(arm, tex)

# ---- LaTeX ---------------------------------------------------------------
latex_table <- function(name) {
  t <- tables[[name]]
  nc <- length(t$cols)
  grp <- vapply(t$cols, `[[`, "", "group")
  # Group header row with \cmidrule under each group (columns offset by Setting and Method).
  runs <- rle(grp)
  ends <- cumsum(runs$lengths) + 2
  starts <- ends - runs$lengths + 1
  top <- c("", "", unlist(Map(function(g, l) if (nzchar(g)) sprintf("\\multicolumn{%d}{c}{%s}", l, t$groups[[g]]) else rep("", l),
                                runs$values, runs$lengths)))
  rules <- paste(sprintf("\\cmidrule(lr){%d-%d}", starts[nzchar(runs$values)], ends[nzchar(runs$values)]), collapse = "")
  heads <- c("Setting", t$label_head, vapply(t$cols, `[[`, "", "head"))
  body <- vapply(seq_along(t$rows), function(i) {
    arm <- t$rows[i]
    paste0(if (i == 1) setting else "", " & ", row_label(name, arm, TRUE), " & ",
           paste(vapply(t$cols, function(c) cell_text(name, arm, c, TRUE), ""), collapse = " & "), " \\\\")
  }, "")
  c("% Generated by inst/scripts/null_sim/12_paper_tables.R; do not edit.",
    sprintf("\\begin{tabular}{ll%s}", strrep("r", nc)), "\\toprule",
    if (any(nzchar(grp))) c(paste(paste(top, collapse = " & "), "\\\\"), rules),
    paste(paste(heads, collapse = " & "), "\\\\"), "\\midrule", body, "\\bottomrule", "\\end{tabular}")
}
settings_rows <- list(
  c("Setting", setting),
  c("Chromosome", "22"),
  c("LD blocks", "fixed 250 kb, shared by both panels"),
  c("GWAS composition", sprintf("%s %.0f\\%%, %s %.0f\\%% (Stouffer meta-analysis)", pa, 100 * w0, pb, 100 * (1 - w0))),
  c("Replicates", sprintf("%d null, %d non-null", length(unique(pm$rep)), length(unique(nn$rep)))),
  c("Mixing weights", sprintf("$w$ on %s: 0.0 to 1.0 in steps of 0.1", pa)),
  c("PRS-CS", sprintf("$\\phi$ %s, sampler defaults (1,000 iterations, 500 burn-in, thinning 5)", phi_mode)),
  c("Sampler seeds", "1000 + replicate; seed rerun 2000 + replicate"))
latex_settings <- c("% Generated by inst/scripts/null_sim/12_paper_tables.R; do not edit.",
                    "\\begin{tabular}{ll}", "\\toprule", "Item & Value \\\\", "\\midrule",
                    vapply(settings_rows[-1], function(r) sprintf("%s & %s \\\\", r[1], r[2]), ""),
                    "\\bottomrule", "\\end{tabular}")
writeLines(latex_table("null"), file.path(out_dir, "table_null.tex"))
writeLines(latex_table("nonnull"), file.path(out_dir, "table_nonnull.tex"))
writeLines(latex_table("seed_contrast"), file.path(out_dir, "table_seed_contrast.tex"))
writeLines(latex_settings, file.path(out_dir, "table_settings.tex"))

# ---- HTML ----------------------------------------------------------------
tex2html <- function(s) {
  s <- gsub("\\\\%", "%", s)
  s <- gsub("vs\\.\\\\ ", "vs. ", s)
  s <- gsub("\\$r\\^2\\$", "<i>r</i><sup>2</sup>", s)
  s <- gsub("\\$\\\\times 10\\^3\\$", "&times;10<sup>3</sup>", s)
  s <- gsub("\\$\\\\phi\\$", "<i>&phi;</i>", s)
  s <- gsub("\\$w = ([0-9.]+)\\$", "<i>w</i> = \\1", s)
  s <- gsub("\\$w\\$", "<i>w</i>", s)
  s
}
html_table <- function(name) {
  t <- tables[[name]]
  grp <- vapply(t$cols, `[[`, "", "group")
  runs <- rle(grp)
  top <- paste0("<th></th><th></th>", paste(unlist(Map(function(g, l) {
    if (nzchar(g)) sprintf("<th colspan=\"%d\" class=\"grp\">%s</th>", l, tex2html(t$groups[[g]])) else rep("<th></th>", l)
  }, runs$values, runs$lengths)), collapse = ""))
  heads <- paste0(sprintf("<th>Setting</th><th>%s</th>", t$label_head),
                  paste(sprintf("<th>%s</th>", vapply(t$cols, function(c) tex2html(c$head), "")), collapse = ""))
  body <- vapply(seq_along(t$rows), function(i) {
    arm <- t$rows[i]
    sprintf("<tr><td>%s</td><td>%s</td>%s</tr>", if (i == 1) setting else "", row_label(name, arm, FALSE),
            paste(sprintf("<td>%s</td>", vapply(t$cols, function(c) cell_text(name, arm, c, FALSE), "")), collapse = ""))
  }, "")
  c("<table>", "<thead>", if (any(nzchar(grp))) sprintf("<tr class=\"top\">%s</tr>", top),
    sprintf("<tr class=\"%s\">%s</tr>", if (any(nzchar(grp))) "sub" else "top", heads), "</thead>",
    "<tbody>", body, "</tbody>", "</table>")
}
html_settings <- c("<table>", "<thead>", "<tr class=\"top\"><th>Setting</th>",
                   paste(sprintf("<th>%s</th>", vapply(settings_rows[-1], `[[`, "", 1)), collapse = ""), "</tr>", "</thead>",
                   "<tbody>", "<tr>", paste(sprintf("<td>%s</td>", tex2html(vapply(settings_rows, `[[`, "", 2))), collapse = ""),
                   "</tr>", "</tbody>", "</table>")
n_null <- length(unique(pm$rep))
n_nonnull <- length(unique(nn$rep))
common_caption <- sprintf(paste(
  "Chromosome 22, fixed 250 kb LD blocks; the GWAS combines %s and %s discovery samples at %.0f:%.0f.",
  "PRS-CS with automatic &phi;. Means over replicates with Monte Carlo standard errors in parentheses; within a replicate",
  "every method analyses the same simulated GWAS and scores the same target individuals. PRS-MIX at <i>w</i> = 1.0 and",
  "0.0 is identical to %s-only and %s-only PRS-CS. The seed rerun repeats %s-only PRS-CS with a different sampler seed",
  "and measures run-to-run variation. <sup>*</sup>Generating weight."), pa, pb, 100 * w0, 100 * (1 - w0), pa, pb, pa)
captions <- list(
  null = sprintf(paste("Null GWAS (no true effects, %d replicates). LD-fit loss is the excess Gaussian loss (ridge 0.001) of each",
                       "method's LD matrix over the mixed matrix at the generating weight; lower is better. Null score variance",
                       "compares the variance of target-sample scores with %s-only PRS-CS for the same GWAS, as a percentage change;",
                       "negative values mean less spread in scores that carry no true signal. Ratios are averaged on the log scale",
                       "(geometric mean) with delta-method standard errors. %s"), n_null, pa, common_caption),
  nonnull = sprintf(paste("Non-null GWAS (sparse causal effects shared across ancestries, %d replicates). Top-10%% overlap is the",
                          "share of target individuals in a method's top score decile who are also in %s-only PRS-CS's top decile;",
                          "correlation is between their scores. <i>r</i><sup>2</sup> is the squared correlation between scores and",
                          "true genetic values. %s"), n_nonnull, pa, common_caption),
  seed_contrast = sprintf(paste("Replicate-paired comparisons of each method with the %s-only seed rerun, with standard errors in parentheses.",
                                "Null score variance is the method's variance relative to the seed rerun's, as a percentage change",
                                "(geometric mean over replicates; negative means less noise than a reseeded %s-only run). Top-10%% overlap",
                                "and correlation are first measured against %s-only PRS-CS, as in Table 2, and the columns give the method",
                                "minus the seed rerun; negative values mean a method moves scores further than a new sampler seed alone",
                                "does. pp, percentage points. <sup>*</sup>Generating weight."), pa, pa, pa),
  settings = "Simulation and analysis settings for this study.")
titles <- c(null = "Table 1. Null GWAS simulation", nonnull = "Table 2. Non-null GWAS simulation",
            seed_contrast = "Table 3. Replicate-paired differences against the seed rerun", settings = "Table 4. Simulation settings")
section <- function(name, tab) c("<section>", sprintf("<h2>%s</h2>", titles[[name]]),
                                 sprintf("<p class=\"cap\">%s</p>", captions[[name]]), "<div class=\"wrap\">", tab, "</div>", "</section>")
head_commit <- git("rev-parse", "HEAD")
html <- c(
  "<!doctype html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
  sprintf("<title>%s tables</title>", basename(normalizePath(out_dir))),
  "<!-- Generated by inst/scripts/null_sim/12_paper_tables.R; do not edit. -->",
  "<style>",
  "body { margin: 0; padding: 40px 24px 56px; background: #fff; color: #111; font: 15px/1.45 Georgia, \"Times New Roman\", Times, serif; }",
  "main { max-width: 1180px; margin: 0 auto; display: flex; flex-direction: column; gap: 44px; }",
  "h2 { font: 600 17px/1.3 -apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif; margin: 0 0 8px; }",
  ".cap { margin: 0 0 14px; font-size: 14px; max-width: none; }",
  ".wrap { overflow-x: auto; }",
  "table { border-collapse: collapse; width: 100%; font-variant-numeric: lining-nums tabular-nums; font-size: 14.5px; }",
  "th, td { padding: 4px 10px; text-align: center; white-space: nowrap; vertical-align: baseline; }",
  "th:nth-child(-n+2), td:nth-child(-n+2) { text-align: left; }",
  "thead th { font-weight: 700; }",
  "thead tr.top th { border-top: 2px solid #111; padding-top: 6px; }",
  "thead tr:last-child th { border-bottom: 1px solid #111; padding-bottom: 6px; }",
  "thead th.grp { font-weight: 700; border-bottom: 1px solid #111; }",
  "tbody tr:last-child td { border-bottom: 2px solid #111; padding-bottom: 6px; }",
  "tbody tr:first-child td { padding-top: 6px; }",
  "footer { font: 12px/1.5 -apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif; color: #666; }",
  "code { font-size: 11.5px; }",
  "</style>", "</head>", "<body>", "<main>",
  section("null", html_table("null")),
  section("nonnull", html_table("nonnull")),
  section("seed_contrast", html_table("seed_contrast")),
  section("settings", html_settings),
  sprintf(paste("<footer>Rendered from <code>%s/inputs/</code> by <code>%s</code> at PRS-MIX %s. The <code>.tex</code> files in",
                "the same folder hold the same values, and <code>tables.tsv</code> holds them unrounded with replicate",
                "counts.</footer>"), out_dir, program, substr(head_commit, 1, 7)),
  "</main>", "</body>", "</html>")
writeLines(html, file.path(out_dir, "tables.html"))

# ---- Manifest ------------------------------------------------------------
out_files <- c(table_null = "table_null.tex", table_nonnull = "table_nonnull.tex",
               table_seed_contrast = "table_seed_contrast.tex", table_settings = "table_settings.tex",
               tables_tsv = "tables.tsv", tables_html = "tables.html")
outputs <- lapply(out_files, function(f) list(path = file.path(out_dir, f), sha256 = sha(file.path(out_dir, f))))
prog_dirty <- length(git("status", "--porcelain", "--", program, "inst/scripts/null_sim/common.R")) > 0
settings <- list(
  chromosome = 22, populations = pops, block_layout = "fixed 250 kb",
  generating_composition = stats::setNames(as.list(c(w0, 1 - w0)), pops),
  mixing_weights = seq(0, 1, by = 0.1), phi_mode = phi_mode,
  replicates = list(null = n_null, nonnull = n_nonnull),
  seed_control = "single-panel arm rerun with sampler seed 2000 + rep (main runs use 1000 + rep)",
  uncertainty = "mean over replicates with Monte Carlo SE (SD / sqrt(n)); paired differences use the SE of the per-replicate difference; variance ratios use the log scale with delta-method SE",
  metric_definitions = list(
    ld_fit_excess = "Gaussian LD-fit loss (ridge delta 0.001) of the method's LD matrix to the generating covariance of the null z-scores, minus the loss at the generating weight, per replicate",
    null_var_pct = sprintf("100 * (exp(mean over replicates of log(variance of target scores / variance under %s-only PRS-CS)) - 1), null GWAS; delta-method SE", pa),
    top10_pct = sprintf("share of target individuals in the method's top score decile also in %s-only PRS-CS's top decile, non-null GWAS", pa),
    cor = sprintf("correlation of target scores with %s-only PRS-CS scores, non-null GWAS", pa),
    r2 = "squared correlation between target scores and true genetic values, non-null GWAS",
    seed_contrast = sprintf("null score variance of the method relative to the %s-only seed rerun (geometric mean percentage change); per-replicate difference between the method and the seed rerun in top10_pct (percentage points) and cor (x1000)", pa)))
manifest <- list(
  study = basename(normalizePath(out_dir)), program = program,
  command = paste("Rscript", program, paste(shQuote(args), collapse = " ")),
  asset_commit = head_commit, program_uncommitted = prog_dirty,
  experiment_commit = unique(vapply(inputs, `[[`, "", "source_commit")),
  r_version = R.version.string,
  inputs = unname(lapply(inputs, function(x) x[c("path", "source", "source_commit", "sha256")])),
  outputs = outputs, settings = settings,
  notes = paste("Replicates missing an arm are dropped from that arm's paired comparisons; n per cell is in tables.tsv.",
                "Experiment runs were made on RedHawk with the null-simulation job scripts (cluster/jobs/null_sim, not version-controlled)."))
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(sprintf("Wrote %s (%d cells) at PRS-MIX %s%s", out_dir, nrow(cells), substr(head_commit, 1, 7),
                if (prog_dirty) " with uncommitted program changes" else ""))
