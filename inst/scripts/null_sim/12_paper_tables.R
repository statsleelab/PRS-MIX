#!/usr/bin/env Rscript
# Curated manuscript tables for one null/non-null simulation study (plan
# section 7), written into paper/<study>/ for import by the manuscript
# repository.
#
# Rows: EUR-only PRS-CS (reference), EUR-only PRS-CS rerun with another sampler
# seed (run-to-run control), PRS-MIX at the given weights, and ASN-only PRS-CS.
# Weights are on pops[1].
#
#   table_null.tex     Null GWAS: excess LD-fit loss over the generating weight
#                      and null score variance relative to EUR-only PRS-CS, in
#                      each target ancestry.
#   table_nonnull.tex  Non-null GWAS: top-decile overlap with EUR-only PRS-CS
#                      and r^2 with the true genetic value, in each target
#                      ancestry.
#   tables.tsv         Every number in both tables, unrounded.
#   inputs/            Frozen copies of the metric tables used.
#   manifest.json      Program, command, commits, input and output sha256.
#
# Each cell is the mean over replicates with a percentile bootstrap 95% CI
# (replicates resampled). Comparisons with EUR-only PRS-CS are paired within
# replicate (same simulated GWAS and target individuals).
#
# Usage (from the PRS-MIX root, with the inputs committed):
#   Rscript inst/scripts/null_sim/12_paper_tables.R <metrics_dir> <out_dir> \
#     [phi_mode=auto] [weights=0.9,0.7,0.5] [w0=0.7] [pops=EUR,ASN] [n_boot=10000]

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
usage <- paste("Rscript 12_paper_tables.R <metrics_dir> <out_dir> [phi_mode=auto]",
               "[weights=0.9,0.7,0.5] [w0=0.7] [pops=EUR,ASN] [n_boot=10000]")
args <- ns_args(usage, 2)
metrics_dir <- args[[1]]
out_dir <- args[[2]]
phi_mode <- ns_arg(args, 3, "auto")
weights <- as.numeric(strsplit(ns_arg(args, 4, "0.9,0.7,0.5"), ",")[[1]])
w0 <- as.numeric(ns_arg(args, 5, "0.7"))
pops <- toupper(strsplit(ns_arg(args, 6, "EUR,ASN"), ",")[[1]])
n_boot <- as.integer(ns_arg(args, 7, "10000"))
boot_seed <- 20261105L
stopifnot(length(pops) == 2)
pa <- pops[1]
pb <- pops[2]

sha <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("\\s.*$", "", out)
}
git <- function(...) system2("git", c(...), stdout = TRUE)

# ---- Freeze inputs -------------------------------------------------------
input_names <- c("prscs_metrics.tsv", "nonnull_metrics.tsv", "ld_fit.tsv")
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
inputs <- lapply(input_names, function(f) {
  src <- file.path(metrics_dir, f)
  dst <- file.path(out_dir, "inputs", f)
  file.copy(src, dst, overwrite = TRUE)
  stopifnot(sha(src) == sha(dst))
  dirty <- length(git("status", "--porcelain", "--", src)) > 0
  list(path = dst, source = src, sha256 = sha(dst),
       source_commit = git("log", "-1", "--format=%H", "--", src),
       source_uncommitted = dirty)
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

# ---- Statistics ----------------------------------------------------------
set.seed(boot_seed)
summ <- function(x) {
  x <- x[is.finite(x)]
  b <- vapply(seq_len(n_boot), function(i) mean(sample(x, replace = TRUE)), numeric(1))
  c(mean = mean(x), ci_lo = unname(stats::quantile(b, 0.025)),
    ci_hi = unname(stats::quantile(b, 0.975)), n = length(x))
}
# Per-replicate value of `col` for `arm`, optionally paired with `ref_arm`.
by_rep <- function(tab, arm, col, ref_arm = NULL, f = function(a, r) a) {
  a <- tab[tab$arm == arm, c("rep", col)]
  if (is.null(ref_arm)) return(a[[2]])
  r <- tab[tab$arm == ref_arm, c("rep", col)]
  z <- merge(a, r, by = "rep")
  f(z[[2]], z[[3]])
}
ref <- paste0("single_", pa)
rows <- c(ref, paste0(ref, "_seed2"), sprintf("mix_w%.1f", weights), paste0("single_", pb))
row_weight <- function(arm) ns_arm_weight(arm, pops)

stats <- list()
add <- function(table, arm, column, target, x, note = "") {
  s <- if (length(x)) summ(x) else c(mean = NA, ci_lo = NA, ci_hi = NA, n = 0)
  stats[[length(stats) + 1]] <<- data.frame(table = table, arm = arm, column = column, target = target,
                                            mean = s[["mean"]], ci_lo = s[["ci_lo"]], ci_hi = s[["ci_hi"]],
                                            n_rep = s[["n"]], note = note)
}
ld0 <- ld[abs(ld$delta - 0.001) < 1e-12, ]
for (arm in rows) {
  seed_row <- grepl("_seed[0-9]+$", arm)
  # Null: excess LD-fit loss L(w) - L(w0), per replicate (does not depend on the sampler).
  if (seed_row) {
    add("null", arm, "ld_fit_excess", "", numeric(0), "not applicable: LD fit does not depend on the sampler seed")
  } else {
    w <- row_weight(arm)
    a <- ld0[abs(ld0$w - w) < 1e-9, c("rep", "L")]
    b <- ld0[abs(ld0$w - w0) < 1e-9, c("rep", "L")]
    z <- merge(a, b, by = "rep")
    add("null", arm, "ld_fit_excess", "", z[[2]] - z[[3]])
  }
  for (p in pops) {
    col <- paste0("score_var_", p)
    if (arm == ref) {
      add("null", arm, "null_score_var_pct", p, numeric(0), "reference")
    } else {
      add("null", arm, "null_score_var_pct", p, by_rep(pm, arm, col, ref, function(a, r) 100 * (a / r - 1)))
    }
    if (arm == ref) {
      add("nonnull", arm, "top10_overlap_pct", p, numeric(0), "reference")
    } else {
      add("nonnull", arm, "top10_overlap_pct", p, 100 * by_rep(nn, arm, sprintf("top10_%s_%s", pa, p)))
    }
    add("nonnull", arm, "r2", p, by_rep(nn, arm, paste0("r2_", p)))
  }
}
stats <- do.call(rbind, stats)

# ---- LaTeX ---------------------------------------------------------------
label <- function(arm) {
  if (arm == ref) return(sprintf("PRS-CS, %s panel", pa))
  if (grepl("_seed[0-9]+$", arm)) return(sprintf("\\quad rerun, new sampler seed"))
  if (grepl("^single_", arm)) return(sprintf("PRS-CS, %s panel", sub("^single_", "", arm)))
  w <- row_weight(arm)
  sprintf("PRS-MIX, $w_{\\mathrm{%s}} = %.1f$%s", pa, w, if (abs(w - w0) < 1e-9) "$^{*}$" else "")
}
num <- function(x, d) {
  s <- formatC(x, format = "f", digits = d)
  ifelse(x < 0, paste0("$-$", sub("^-", "", s)), s)
}
cell <- function(table, arm, column, target, d) {
  r <- stats[stats$table == table & stats$arm == arm & stats$column == column & stats$target == target, ]
  if (r$note == "reference") return("ref.")
  if (is.na(r$mean)) return("--")
  if (r$mean == 0 && r$ci_lo == 0 && r$ci_hi == 0) return("0")  # zero by definition (generating weight)
  sprintf("%s {\\scriptsize[%s, %s]}", num(r$mean, d), num(r$ci_lo, d), num(r$ci_hi, d))
}
tabular <- function(table, cols, header) {
  body <- vapply(rows, function(arm) {
    paste0(label(arm), " & ", paste(vapply(cols, function(cc) cell(table, arm, cc[[1]], cc[[2]], cc[[3]]), ""), collapse = " & "), " \\\\")
  }, "")
  c(sprintf("%% Generated by inst/scripts/null_sim/12_paper_tables.R; do not edit."),
    sprintf("\\begin{tabular}{l%s}", strrep("r", length(cols))),
    "\\toprule", header, "\\midrule", body, "\\bottomrule", "\\end{tabular}")
}
null_cols <- list(list("ld_fit_excess", "", 3), list("null_score_var_pct", pa, 1), list("null_score_var_pct", pb, 1))
null_tex <- tabular("null", null_cols, c(
  sprintf(" & LD-fit loss & \\multicolumn{2}{c}{Null score variance vs.\\ %s-only (\\%%)} \\\\", pa),
  "\\cmidrule(lr){3-4}",
  sprintf("Method & above $w = %.1f$ & %s targets & %s targets \\\\", w0, pa, pb)))
nonnull_cols <- list(list("top10_overlap_pct", pa, 1), list("top10_overlap_pct", pb, 1), list("r2", pa, 3), list("r2", pb, 3))
nonnull_tex <- tabular("nonnull", nonnull_cols, c(
  sprintf(" & \\multicolumn{2}{c}{Top-10\\%% overlap with %s-only (\\%%)} & \\multicolumn{2}{c}{$r^2$ with true genetic value} \\\\", pa),
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  sprintf("Method & %s targets & %s targets & %s targets & %s targets \\\\", pa, pb, pa, pb)))
writeLines(null_tex, file.path(out_dir, "table_null.tex"))
writeLines(nonnull_tex, file.path(out_dir, "table_nonnull.tex"))
ns_write_tsv(stats, file.path(out_dir, "tables.tsv"))

# ---- Manifest ------------------------------------------------------------
program <- "inst/scripts/null_sim/12_paper_tables.R"
head <- git("rev-parse", "HEAD")
prog_dirty <- length(git("status", "--porcelain", "--", program, "inst/scripts/null_sim/common.R")) > 0
cmd <- paste("Rscript", program, paste(shQuote(args), collapse = " "))
outputs <- lapply(c("table_null.tex", "table_nonnull.tex", "tables.tsv"), function(f) {
  list(path = file.path(out_dir, f), sha256 = sha(file.path(out_dir, f)))
})
names(outputs) <- c("table_null", "table_nonnull", "tables_tsv")
settings <- list(
  chromosome = 22, populations = pops, block_layout = "fixed 250 kb",
  generating_composition = stats::setNames(as.list(c(w0, 1 - w0)), pops),
  mixing_weights = weights, phi_mode = phi_mode,
  replicates = list(null = length(unique(pm$rep)), nonnull = length(unique(nn$rep))),
  seed_control = "single-panel arm rerun with sampler seed 2000 + rep (main runs use 1000 + rep)",
  uncertainty = sprintf("mean over replicates, percentile bootstrap 95%% CI (%d resamples of replicates, seed %d)", n_boot, boot_seed),
  metric_definitions = list(
    ld_fit_excess = "Gaussian LD-fit loss (ridge delta 0.001) of the arm's LD matrix to the generating covariance of the null z-scores, minus the loss at the generating weight; per replicate, then averaged",
    null_score_var_pct = sprintf("100 * (variance of target scores / variance under %s-only PRS-CS - 1), null GWAS, paired by replicate", pa),
    top10_overlap_pct = sprintf("share of target individuals in the arm's top score decile who are also in %s-only PRS-CS's top decile, non-null GWAS", pa),
    r2 = "squared correlation between target scores and true genetic values, non-null GWAS"
  )
)
manifest <- list(
  study = basename(normalizePath(out_dir)),
  program = program, command = cmd,
  asset_commit = head, program_uncommitted = prog_dirty,
  experiment_commit = unique(vapply(inputs, `[[`, "", "source_commit")),
  r_version = R.version.string,
  inputs = unname(lapply(inputs, function(x) x[c("path", "source", "source_commit", "sha256")])),
  outputs = outputs, settings = settings,
  notes = paste("Replicates missing an arm are dropped from that arm's paired comparisons; n_rep per cell is in tables.tsv.",
                "Experiment runs were made on RedHawk with the null-simulation job scripts (cluster/jobs/null_sim, not version-controlled).")
)
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(sprintf("Wrote %s (%d table cells) at PRS-MIX %s%s", out_dir, nrow(stats), substr(head, 1, 7),
                if (prog_dirty) " with uncommitted program changes" else ""))
