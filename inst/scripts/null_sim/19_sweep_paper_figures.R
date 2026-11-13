#!/usr/bin/env Rscript
# Paper figure for the 250 kb ancestry-pair x GWAS-composition sweep: who is in
# the top 10%, PRS-MIX against the EUR panel's own rerun.
#
#   fig_top10_overlap_eur.{png,pdf}
#     For the three EUR pairs (EUR+ASN, EUR+AFR, EUR+SAS), scored in EUR
#     (top row) and in the other ancestry (bottom row): per-replicate overlap
#     with the EUR panel's top 10% for
#       grey  the EUR panel rerun with a new sampler seed (2000 + rep), the
#             agreement expected from sampler randomness alone;
#       blue  PRS-MIX at the generating composition (seed 1000 + rep, like
#             the EUR panel run it is compared with).
#     Boxes over replicates with every replicate drawn as a point.
#   figure_data.tsv   every plotted value (one row per replicate x box)
#   inputs/           frozen copies of the metric tables used
#   manifest.json     program, command, commits, input and output sha256
#
# ASN+AFR is left out so that every panel has the same reference (the EUR
# panel); its top-10% results are in 18_sweep_paper_tables.R's Table 2.
#
# Usage (from the PRS-MIX root):
#   Rscript inst/scripts/null_sim/19_sweep_paper_figures.R <sweep_dir> <out_dir> [reps=all] [phi_mode=auto] [provisional=no]
#
#   <sweep_dir>  e.g. inst/scripts/results/sim250_sweep_100rep, holding
#                <pair>/w<w>/metrics/nonnull_metrics.tsv for every setting.
#   reps         all, or N to use replicates 1..N of every setting.
#
# As in 18_sweep_paper_tables.R, writing under paper/ requires committed inputs
# and program and at least 100 replicates per setting, except for a placeholder:
# provisional=yes and an output path containing "provisional" (e.g.
# paper/<study>_provisional/figures).

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages(library(ggplot2))
grDevices::pdf(NULL)  # keep ggplot from opening a default device (Rplots.pdf) in the working directory
args <- ns_args("Rscript 19_sweep_paper_figures.R <sweep_dir> <out_dir> [reps=all] [phi_mode=auto] [provisional=no]", 2)
sweep_dir <- args[[1]]
out_dir <- args[[2]]
reps_arg <- ns_arg(args, 3, "all")
phi_mode <- ns_arg(args, 4, "auto")
provisional <- ns_arg(args, 5, "no") == "yes"
program <- "inst/scripts/null_sim/19_sweep_paper_figures.R"
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
to_paper <- grepl("(^|/)paper/", normalizePath(out_dir))

sha <- function(path) sub("\\s.*$", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
git <- function(...) system2("git", c(...), stdout = TRUE)

pairs <- list(eur_asn = c("EUR", "ASN"), eur_afr = c("EUR", "AFR"), eur_sas = c("EUR", "SAS"))
ws <- c(0.9, 0.7, 0.5)
settings <- do.call(rbind, lapply(names(pairs), function(p) data.frame(pair = p, w = ws)))

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

# ---- Figure data -----------------------------------------------------------
m_rerun <- "EUR panel alone, rerun with a new seed"
m_mix <- "PRS-MIX (blended panel)"
rows <- list()
n_reps <- integer(0)
for (i in seq_len(nrow(settings))) {
  s <- settings[i, ]
  A <- pairs[[s$pair]][1]; B <- pairs[[s$pair]][2]
  nn <- ns_read_tsv(inputs[[i]]$path)
  nn <- nn[nn$phi_mode == phi_mode, ]
  if (reps_arg != "all") nn <- nn[nn$rep <= as.integer(reps_arg), ]
  mix <- sprintf("mix_w%.1f", s$w)
  rerun <- paste0("single_", A, "_seed2")
  reps <- sort(Reduce(intersect, lapply(c(mix, rerun), function(a) nn$rep[nn$arm == a])))
  label <- sprintf("%s + %s, %.0f:%.0f", A, B, 100 * s$w, 100 * (1 - s$w))
  n_reps[[label]] <- length(reps)
  v <- function(arm, col) { y <- nn[nn$arm == arm, ]; y[[col]][match(reps, y$rep)] }
  for (k in 1:2) {
    tg <- c(A, B)[k]
    col <- sprintf("top10_%s_%s", A, tg)
    for (m in list(c(rerun, m_rerun), c(mix, m_mix))) {
      rows[[length(rows) + 1]] <- data.frame(
        pair = paste(A, "and", B), w = s$w, mix = sprintf("%.0f:%.0f", 100 * s$w, 100 * (1 - s$w)),
        target = tg, row = k, arm = m[1], method = m[2], rep = reps, overlap_pct = 100 * v(m[1], col))
    }
  }
}
fd <- do.call(rbind, rows)
placeholder <- provisional && grepl("provisional", normalizePath(out_dir))
if (provisional && !placeholder) stop("provisional=yes needs an output path containing 'provisional'.", call. = FALSE)
if (to_paper && any(n_reps < 100) && !placeholder) {
  stop(sprintf("Writing under paper/ needs >= 100 replicates per setting; fewest is %d (%s).",
               min(n_reps), names(which.min(n_reps))), call. = FALSE)
}
ns_write_tsv(fd, file.path(out_dir, "figure_data.tsv"))
rep_text <- if (length(unique(n_reps)) == 1) as.character(n_reps[[1]]) else sprintf("%d to %d", min(n_reps), max(n_reps))

# ---- Figure ------------------------------------------------------------------
ink <- "#1f2328"; muted <- "#5f6670"; grid <- "#e6e8eb"
fd$mix <- factor(fd$mix, c("90:10", "70:30", "50:50"))
fd$method <- factor(fd$method, c(m_rerun, m_mix))
fd$strip <- paste0(fd$pair, ", scored in ", fd$target)
fd$strip <- factor(fd$strip, unique(fd$strip[order(fd$row, factor(fd$pair, paste("EUR and", c("ASN", "AFR", "SAS"))))]))
prov_note <- if (any(n_reps < 100)) sprintf("Provisional: %s replicates per setting.", rep_text) else NULL
fig <- ggplot(fd, aes(mix, overlap_pct, fill = method)) +
  geom_boxplot(colour = ink, width = 0.66, linewidth = 0.3, outlier.shape = NA, position = position_dodge(width = 0.76)) +
  geom_point(aes(group = method), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.76, seed = 1),
             size = 0.5, colour = ink, alpha = 0.5) +
  facet_wrap(~strip, nrow = 2) +
  scale_fill_manual(values = setNames(c("#d9dde2", "#6da7ec"), c(m_rerun, m_mix)), name = NULL) +
  scale_y_continuous(limits = c(50, 100), breaks = seq(50, 100, 10)) +
  labs(x = "GWAS ancestry mix (as named in the panel title, e.g. 70:30 = 70% EUR, 30% ASN)",
       y = "Overlap with the EUR panel's top 10% (%)",
       caption = paste(c(
         sprintf("Overlap: share of a run's top 10%% that is also in the EUR panel's original top 10%%. Each point is one replicate (%s per box).", rep_text),
         "Grey: the EUR panel rerun with a new seed, the agreement expected from sampler randomness alone.",
         "Blue below grey: PRS-MIX changes the top 10% more than randomness does.", prov_note), collapse = "\n")) +
  theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = 0.3, colour = ink), axis.ticks = element_line(linewidth = 0.3, colour = ink),
        axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = grid),
        strip.background = element_blank(), strip.text = element_text(face = "bold", colour = ink, hjust = 0),
        legend.position = "top", legend.justification = "left", legend.key.width = unit(14, "pt"),
        plot.caption = element_text(hjust = 0, colour = muted, size = 7))
stem <- file.path(out_dir, "fig_top10_overlap_eur")
ggsave(paste0(stem, ".png"), fig, width = 7.0, height = 5.0, dpi = 300, device = ragg::agg_png)
ggsave(paste0(stem, ".pdf"), fig, width = 7.0, height = 5.0, device = grDevices::pdf, useDingbats = FALSE)

# ---- Manifest ------------------------------------------------------------
head_commit <- git("rev-parse", "HEAD")
out_files <- c(fig_png = "fig_top10_overlap_eur.png", fig_pdf = "fig_top10_overlap_eur.pdf", figure_data = "figure_data.tsv")
outputs <- lapply(out_files, function(f) list(path = file.path(out_dir, f), sha256 = sha(file.path(out_dir, f))))
manifest <- list(
  study = basename(normalizePath(out_dir)), program = program,
  command = paste("Rscript", program, paste(shQuote(args), collapse = " ")),
  asset_commit = head_commit, program_uncommitted = prog_dirty, inputs_uncommitted = inputs_dirty,
  provisional = any(n_reps < 100),
  experiment_commit = unique(vapply(inputs, `[[`, "", "source_commit")),
  r_version = R.version.string, ggplot2_version = as.character(utils::packageVersion("ggplot2")),
  inputs = unname(lapply(inputs, function(x) x[c("path", "source", "source_commit", "sha256")])),
  outputs = outputs,
  settings = list(
    chromosome = 22, block_layout = "fixed 250 kb", pairs = unname(lapply(pairs, as.list)), compositions_A = ws,
    replicates_used = as.list(n_reps), replicate_filter = reps_arg, phi_mode = phi_mode,
    overlap = "top10_EUR_<target> of 07_metrics_nonnull.R: share of a run's top 10% of target scores also in the EUR panel run's top 10%",
    rerun = "single_EUR_seed2: the EUR panel with sampler seed 2000 + rep (all other arms use 1000 + rep)",
    prs_mix = "mix_w<w> at the generating composition w"),
  notes = "Experiment runs were made on RedHawk with the null-simulation job scripts (cluster/jobs/null_sim, not version-controlled).")
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(sprintf("Wrote %s (%d plotted values, %s replicates) at PRS-MIX %s%s", out_dir, nrow(fd), rep_text, substr(head_commit, 1, 7),
                if (inputs_dirty || prog_dirty) " with uncommitted program or inputs" else ""))
