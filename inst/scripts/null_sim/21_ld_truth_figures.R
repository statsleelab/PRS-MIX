#!/usr/bin/env Rscript
# Figures: the LD the null GWAS was simulated under, against the LD each method
# uses, for one 250 kb block, at every GWAS composition. Draws from the files
# written on RedHawk by 20_ld_truth_extract.R.
#
#   fig_ld_truth_<pair>_block<b>.{png,pdf}
#     Rows: GWAS compositions 90:10, 70:30, 50:50 (A : B). Columns:
#       true GWAS LD      w R_A,disc + (1 - w) R_B,disc (Cov(z) of the null z-scores)
#       PRS-MIX LD        w R_A,panel + (1 - w) R_B,panel
#       A panel LD        R_A,panel (standard PRS-CS with the majority-ancestry panel)
#       true - PRS-MIX, true - A panel
#     One signed colour scale (-1 to 1) for every heatmap; lower triangle only.
#   figure_data.tsv.gz    every plotted value (pair, w, column, i, j, value), gzipped
#   block_summary.tsv     per pair and w: mean |true - candidate| over the drawn
#                         block's off-diagonal pairs, and the block's rank among
#                         all chromosome blocks by the A-panel misfit
#   inputs/               frozen copies of the 20_ld_truth_extract.R outputs used
#   manifest.json         program, command, commits, input and output sha256
#
# Usage (from the PRS-MIX root):
#   Rscript inst/scripts/null_sim/21_ld_truth_figures.R <extract_dir> <out_dir> <block> [pairs=eur_asn,eur_afr,eur_sas]
#
#   <extract_dir> may be <out_dir>/inputs itself (then the files are not copied).
#
# These figures use no PRS-CS output and no replicate: the true LD comes from
# the discovery genotypes and the candidates from the reference panels.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages(library(ggplot2))
grDevices::pdf(NULL)  # keep ggplot from opening a default device (Rplots.pdf) in the working directory
args <- ns_args("Rscript 21_ld_truth_figures.R <extract_dir> <out_dir> <block> [pairs=eur_asn,eur_afr,eur_sas]", 3)
extract_dir <- args[[1]]
out_dir <- args[[2]]
block <- as.integer(args[[3]])
pairs <- strsplit(ns_arg(args, 4, "eur_asn,eur_afr,eur_sas"), ",")[[1]]
program <- "inst/scripts/null_sim/21_ld_truth_figures.R"
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
to_paper <- grepl("(^|/)paper/", normalizePath(out_dir))

sha <- function(path) sub("\\s.*$", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
git <- function(...) system2("git", c(...), stdout = TRUE)

# ---- Freeze inputs -------------------------------------------------------
files <-unlist(lapply(pairs, function(p) c(sprintf("%s_block%d_ld.tsv", p, block), sprintf("%s_block_scan.tsv", p), sprintf("%s_params.tsv", p))))
inputs <- lapply(files, function(f) {
  src <- file.path(extract_dir, f)
  dst <- file.path(out_dir, "inputs", f)
  if (normalizePath(src) != normalizePath(dst, mustWork = FALSE)) file.copy(src, dst, overwrite = TRUE)
  stopifnot(sha(src) == sha(dst))
  list(path = dst, source = src, sha256 = sha(dst),
       source_commit = c(git("log", "-1", "--format=%H", "--", src), NA_character_)[1],
       source_uncommitted = length(git("status", "--porcelain", "--", src)) > 0)
})
names(inputs) <- files
inputs_dirty <- any(vapply(inputs, `[[`, TRUE, "source_uncommitted"))
prog_dirty <- length(git("status", "--porcelain", "--", program, "inst/scripts/null_sim/common.R",
                         "inst/scripts/null_sim/20_ld_truth_extract.R")) > 0
if (to_paper && (inputs_dirty || prog_dirty)) {
  stop("Writing under paper/ needs committed inputs and programs; commit them first.", call. = FALSE)
}
inp <- function(f) ns_read_tsv(inputs[[f]]$path)

# ---- Figure data -----------------------------------------------------------
ws <- c(0.9, 0.7, 0.5)
full <- function(d, col, n) { m <- matrix(0, n, n); m[cbind(d$i, d$j)] <- d[[col]]; m[cbind(d$j, d$i)] <- d[[col]]; m }
to_long <- function(m) { n <- nrow(m); k <- which(row(m) <= col(m)); data.frame(i = row(m)[k], j = col(m)[k], value = m[k]) }
off <- function(m) mean(abs(m[upper.tri(m)]))
fd <- list(); summ <- list(); meta <- list()
for (p in pairs) {
  pops <- toupper(strsplit(p, "_")[[1]])
  d <- inp(sprintf("%s_block%d_ld.tsv", p, block))
  n <- max(d$j)
  da <- full(d, "r_disc_a", n); db <- full(d, "r_disc_b", n)
  ra <- full(d, "r_panel_a", n); rb <- full(d, "r_panel_b", n)
  scan <- inp(sprintf("%s_block_scan.tsv", p))
  prm <- inp(sprintf("%s_params.tsv", p))
  cols <- c("True GWAS LD", "PRS-MIX LD", sprintf("%s panel LD", pops[1]), "True - PRS-MIX", sprintf("True - %s panel", pops[1]))
  for (w in ws) {
    truth <- w * da + (1 - w) * db
    mix <- w * ra + (1 - w) * rb
    m <- list(truth, mix, ra, truth - mix, truth - ra)
    for (k in seq_along(m)) fd[[length(fd) + 1]] <- cbind(pair = p, w = w, column = cols[k], to_long(m[[k]]))
    s <- scan[abs(scan$w - w) < 1e-9, ]
    summ[[length(summ) + 1]] <- data.frame(
      pair = paste(pops, collapse = "+"), w = w, block = block, n_snps = n,
      bp_start = min(d$bp_i), bp_end = max(d$bp_j),
      mad_prsmix = off(truth - mix), mad_a_only = off(truth - ra), mad_b_only = off(truth - rb),
      rank_a_only = rank(-s$mad_a_only)[s$block == block], n_blocks = nrow(s))
  }
  meta[[p]] <- list(pops = pops, cols = cols, n = n,
                    n_disc = stats::setNames(as.integer(prm$value[grepl("^n_discovery_", prm$param)]), pops))
}
fd <- do.call(rbind, fd)
summ <- do.call(rbind, summ)
# Every plotted value; gzipped and rounded to 6 significant digits (one row per heatmap cell, ~430k rows).
fd_out <- fd; fd_out$value <- signif(fd_out$value, 6)
gz <- gzfile(file.path(out_dir, "figure_data.tsv.gz"), "w")
utils::write.table(fd_out, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
ns_write_tsv(summ, file.path(out_dir, "block_summary.tsv"))

# ---- Figures ---------------------------------------------------------------
ink <- "#1f2328"; muted <- "#5f6670"
out_files <- c(figure_data = "figure_data.tsv.gz", block_summary = "block_summary.tsv")
for (p in pairs) {
  mt <- meta[[p]]; A <- mt$pops[1]; B <- mt$pops[2]
  x <- fd[fd$pair == p, ]
  x$column <- factor(x$column, mt$cols)
  x$mix <- factor(sprintf("%.0f:%.0f\n(%s : %s)", 100 * x$w, 100 * (1 - x$w), A, B),
                  sprintf("%.0f:%.0f\n(%s : %s)", 100 * ws, 100 * (1 - ws), A, B))
  s <- summ[summ$pair == paste(A, B, sep = "+"), ]
  fig <- ggplot(x, aes(i, j, fill = value)) + geom_raster() +
    facet_grid(mix ~ column, switch = "y") +
    scale_fill_gradient2(low = "#104281", mid = "#f0efec", high = "#a3201b", midpoint = 0, limits = c(-1, 1),
                         breaks = c(-1, -0.5, 0, 0.5, 1), oob = scales::squish, name = expression("LD " * italic(r))) +
    scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + coord_fixed() +
    labs(title = sprintf("%s + %s", A, B), caption = paste0(
      sprintf("Null simulation, chr22, 250 kb blocks: block %d (%.2f-%.2f Mb, %d SNPs). True GWAS LD: w x %s + (1 - w) x %s from the discovery cohorts\n",
              block, s$bp_start[1] / 1e6, s$bp_end[1] / 1e6, mt$n, A, B),
      sprintf("(n %d / %d), the LD the null z-scores were simulated with. PRS-MIX LD: the same blend of the reference panels; %s panel LD: standard PRS-CS.\n",
              mt$n_disc[[A]], mt$n_disc[[B]], A),
      "Mean |difference| from the true LD, off the diagonal:\n",
      paste0(sprintf("    %.0f:%.0f   PRS-MIX %.3f;  %s panel alone %.3f;  %s panel alone %.3f   (block rank by %s-panel misfit: %d of %d)",
                     100 * s$w, 100 * (1 - s$w), s$mad_prsmix, A, s$mad_a_only, B, s$mad_b_only, A, s$rank_a_only, s$n_blocks),
             collapse = "\n"))) +
    theme_void(base_size = 9, base_family = "Helvetica") +
    theme(plot.title = element_text(face = "bold", colour = ink, size = 10, margin = margin(b = 6)),
          strip.text.x = element_text(face = "bold", colour = ink, size = 8.5, margin = margin(b = 4)),
          strip.text.y.left = element_text(colour = ink, size = 8, angle = 0, hjust = 1, margin = margin(r = 6)),
          panel.spacing.x = unit(6, "pt"), panel.spacing.y = unit(8, "pt"),
          legend.title = element_text(colour = muted, size = 8), legend.text = element_text(colour = muted, size = 7),
          plot.caption = element_text(hjust = 0, colour = muted, size = 6.5, margin = margin(t = 8)),
          plot.margin = margin(8, 8, 8, 8))
  stem <- sprintf("fig_ld_truth_%s_block%d", p, block)
  ggsave(file.path(out_dir, paste0(stem, ".png")), fig, width = 22, height = 14.5, units = "cm", dpi = 300, bg = "white",
         device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), fig, width = 22, height = 14.5, units = "cm", device = grDevices::pdf,
         useDingbats = FALSE)
  out_files[[paste0(p, "_png")]] <- paste0(stem, ".png")
  out_files[[paste0(p, "_pdf")]] <- paste0(stem, ".pdf")
}

# ---- Manifest ------------------------------------------------------------
head_commit <- git("rev-parse", "HEAD")
outputs <- lapply(out_files, function(f) list(path = file.path(out_dir, f), sha256 = sha(file.path(out_dir, f))))
manifest <- list(
  study = basename(normalizePath(out_dir)), program = program,
  command = paste("Rscript", program, paste(shQuote(args), collapse = " ")),
  asset_commit = head_commit, program_uncommitted = prog_dirty, inputs_uncommitted = inputs_dirty,
  experiment_commit = unique(stats::na.omit(vapply(inputs, `[[`, "", "source_commit"))),
  r_version = R.version.string, ggplot2_version = as.character(utils::packageVersion("ggplot2")),
  inputs = unname(lapply(inputs, function(x) x[c("path", "source", "source_commit", "sha256")])),
  outputs = outputs,
  settings = list(
    chromosome = 22, block_layout = "fixed 250 kb", block = block, pairs = pairs, compositions_A = ws,
    true_ld = "w R_A,disc + (1 - w) R_B,disc, R_k,disc = Z_k'Z_k / n_k on the discovery individuals (04_simulate_null.R)",
    prs_mix = "w R_A,panel + (1 - w) R_B,panel from the reference-set panels",
    summary = "mean |true - candidate| over off-diagonal SNP pairs of the block"),
  notes = paste("Inputs were written on RedHawk by 20_ld_truth_extract.R from the sim250 null-simulation directories",
                "(<pair>/w0.7: split, panels and discovery genotypes shared by every composition of the pair)."))
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(sprintf("Wrote %s (block %d, %s) at PRS-MIX %s%s", out_dir, block, paste(pairs, collapse = ", "), substr(head_commit, 1, 7),
                if (inputs_dirty || prog_dirty) " with uncommitted programs or inputs" else ""))
