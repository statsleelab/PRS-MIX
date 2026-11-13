#!/usr/bin/env Rscript
# Figures from 08_summarize.R's tables.
#
#   fig_weight_sweep.png  A: LD fit relative to the generating weight (null);
#                         B: null score variance and C: non-null r^2 in target
#                         individuals by ancestry, against the EUR weight given
#                         to PRS-MIX. The sweep ends are the single-population
#                         arms (identical to vanilla PRS-CS, checked per run).
#   fig_contrasts.png     Primary contrasts of PRS-MIX at the generating weight
#                         against each single-population arm, with 95% CIs.
#
# Usage:
#   Rscript 09_plot.R <metrics_dir> [out_dir=<metrics_dir>] [w0=0.7] [delta=0.001] [phi=auto]

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

args <- ns_args("Rscript 09_plot.R <metrics_dir> [out_dir] [w0=0.7] [delta=0.001] [phi=auto]", 1)
mdir <- args[[1]]
out_dir <- ns_arg(args, 2, mdir)
w0 <- as.numeric(ns_arg(args, 3, "0.7"))
delta0 <- as.numeric(ns_arg(args, 4, "0.001"))
phi0 <- ns_arg(args, 5, "auto")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

curves <- ns_read_tsv(file.path(mdir, "summary_curves.tsv"))
con <- ns_read_tsv(file.path(mdir, "summary_contrasts.tsv"))
n_rep <- max(con$n_rep)

blue <- "#2a78d6"; orange <- "#eb6834"; ink <- "#0b0b0b"; ink2 <- "#52514e"
muted <- "#898781"; grid <- "#e1e0d9"; surf <- "#fcfcfb"
anc_cols <- c("EUR patients" = blue, "ASN patients" = orange)

base_theme <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = surf, colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = grid, linewidth = 0.3),
        axis.text = element_text(colour = muted), axis.title = element_text(colour = ink2, size = 9.5),
        plot.title = element_text(colour = ink, face = "bold", size = 10.5))

x_scale <- scale_x_continuous(
  breaks = c(0, 0.25, 0.5, w0, 1),
  labels = c("0\nASN only", "0.25", "0.5", sprintf("%g\ntrue mix", w0), "1\nEUR only"),
  expand = expansion(mult = c(0.05, 0.05)))
x_lab <- "EUR weight given to PRS-MIX"
truth_line <- geom_vline(xintercept = w0, linetype = "22", colour = muted, linewidth = 0.5)

# Single series (LD fit: no target ancestry), on a pseudo-log axis so both the
# large single-panel values and the small differences near w0 stay readable.
sweep_panel_ld <- function(d, title, ylab) {
  ggplot(d, aes(w, mean)) + truth_line +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = ink2, alpha = 0.18) +
    geom_line(linewidth = 0.8, colour = ink2) +
    geom_point(size = 2.1, shape = 21, stroke = 0.8, fill = surf, colour = ink2) +
    x_scale + scale_y_continuous(transform = scales::pseudo_log_trans(sigma = 0.01, base = 10),
                                 breaks = c(0, 0.01, 0.1, 1, 4)) +
    labs(title = title, x = x_lab, y = ylab) + base_theme
}

sweep_panel_anc <- function(d, title, ylab) {
  ggplot(d, aes(w, mean, colour = anc, fill = anc)) + truth_line +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), colour = NA, alpha = 0.18) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.1, shape = 21, stroke = 0.8, fill = surf) +
    scale_colour_manual(values = anc_cols, name = NULL) + scale_fill_manual(values = anc_cols, name = NULL) +
    x_scale + labs(title = title, x = x_lab, y = ylab) + base_theme
}

ld <- curves[curves$metric == sprintf("ld_fit_dL_delta%g", delta0), ]
pA <- sweep_panel_ld(ld, "A  LD fit to the null GWAS",
                     sprintf("LD-fit score minus score at w = %g\n(lower = better fit; pseudo-log axis)", w0))

anc_curve <- function(source, prefix, scale = 1) {
  d <- curves[curves$source == source & curves$phi_mode == phi0 & curves$metric %in% paste0(prefix, c("EUR", "ASN")), ]
  d$anc <- factor(ifelse(grepl("EUR$", d$metric), "EUR patients", "ASN patients"), names(anc_cols))
  d[, c("mean", "ci_lo", "ci_hi")] <- d[, c("mean", "ci_lo", "ci_hi")] * scale
  d
}
pB <- sweep_panel_anc(anc_curve("null", "score_var_", 1e4), "B  Null: noise in target scores",
                      "Variance of target individuals' scores (×10⁻⁴)\n(lower = less noise)")
pC <- sweep_panel_anc(anc_curve("nonnull", "r2_"), "C  Non-null: recovery of genetic value",
                      "r² between score and true genetic value\n(higher = better)")

fig1 <- (pA | pB | pC) + plot_layout(guides = "collect") +
  plot_annotation(
    title = "PRS-MIX weight sweep, chr22",
    subtitle = sprintf(paste0("GWAS simulated at a true %g:%g EUR:ASN meta-analysis, N = 200,000; PRS-CS %s; ",
                              "LD fit with δ = %g. Points: mean over %d replicates; bands: 95%% bootstrap CI."),
                       100 * w0, 100 * (1 - w0), if (phi0 == "auto") "auto φ" else paste("φ =", phi0), delta0, n_rep),
    theme = theme(plot.background = element_rect(fill = surf, colour = NA),
                  plot.title = element_text(colour = ink, face = "bold"),
                  plot.subtitle = element_text(colour = ink2, size = 9))) &
  theme(legend.position = "bottom", legend.text = element_text(colour = ink2))
ggsave(file.path(out_dir, "fig_weight_sweep.png"), fig1, width = 13, height = 5, dpi = 150, bg = surf)

# Contrasts: one row per (comparison, phi), one panel per metric on its own scale.
prim <- con[con$family == "primary", ]
prim <- prim[!(grepl("^ld_fit", prim$metric) & prim$phi_mode != phi0), ]  # LD fit does not depend on phi
prim$est <- ifelse(prim$scale == "log_ratio", 100 * (exp(prim$mean) - 1), prim$mean)
prim$lo <- ifelse(prim$scale == "log_ratio", 100 * (exp(prim$ci_lo) - 1), prim$ci_lo)
prim$hi <- ifelse(prim$scale == "log_ratio", 100 * (exp(prim$ci_hi) - 1), prim$ci_hi)
metric_lab <- c(
  ld_fit = sprintf("LD fit (δ = %g)\ndifference, lower = better", delta0),
  null_score_var_EUR = "Null score variance, EUR patients\n% change, lower = better",
  null_score_var_ASN = "Null score variance, ASN patients\n% change, lower = better",
  nonnull_r2_EUR = "Non-null r², EUR patients\ndifference, higher = better",
  nonnull_r2_ASN = "Non-null r², ASN patients\ndifference, higher = better")
prim$key <- ifelse(grepl("^ld_fit", prim$metric), "ld_fit", prim$metric)
prim$metric_lab <- factor(metric_lab[prim$key], metric_lab)
prim$anc <- ifelse(grepl("EUR$", prim$metric), "EUR patients", ifelse(grepl("ASN$", prim$metric), "ASN patients", "All SNPs"))
prim$row <- factor(sprintf("vs %s-only", sub("^single_", "", prim$vs)), c("vs ASN-only", "vs EUR-only"))
prim$phi_lab <- factor(ifelse(prim$phi_mode == "auto", "PRS-CS auto φ", paste("PRS-CS φ =", prim$phi_mode)))
prim$sig <- ifelse(!is.na(prim$p_holm) & prim$p_holm < 0.05, "Holm p < 0.05", "not significant")

fig2 <- ggplot(prim, aes(est, row, colour = anc, shape = phi_lab)) +
  geom_vline(xintercept = 0, colour = muted, linewidth = 0.5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.7,
                position = position_dodge(width = 0.5)) +
  geom_point(aes(alpha = sig), size = 2.6, stroke = 0.9, fill = surf, position = position_dodge(width = 0.5)) +
  facet_wrap(~metric_lab, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = c(anc_cols, "All SNPs" = ink2), breaks = names(anc_cols), name = NULL) +
  scale_shape_manual(values = c(16, 21), name = NULL) +
  scale_alpha_manual(values = c("Holm p < 0.05" = 1, "not significant" = 0.45), guide = "none") +
  labs(x = sprintf("PRS-MIX at w = %g minus single-population PRS-CS (95%% bootstrap CI)", w0), y = NULL,
       title = sprintf("PRS-MIX at the true weight vs single-population LD, chr22, %d replicates", n_rep),
       subtitle = paste("Within-replicate paired differences. Faded: not significant after Holm correction.",
                        "LD fit (grey) does not involve PRS-CS.", sep = "\n")) +
  base_theme +
  theme(panel.grid.major.x = element_line(colour = grid, linewidth = 0.3), panel.grid.major.y = element_blank(),
        strip.text = element_text(colour = ink, hjust = 0, size = 9.5),
        plot.subtitle = element_text(colour = ink2, size = 9), legend.position = "bottom",
        legend.box = "vertical", legend.text = element_text(colour = ink2))
ggsave(file.path(out_dir, "fig_contrasts.png"), fig2, width = 8, height = 10, dpi = 150, bg = surf)
message(sprintf("Wrote fig_weight_sweep.png and fig_contrasts.png to %s", out_dir))
