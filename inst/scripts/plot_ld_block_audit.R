#!/usr/bin/env Rscript
# Figures for the LD block-boundary audit, drawn from the CSV tables written by
# summarize_ld_block_audit.R (so they can be made on any machine with a PNG
# device; RedHawk's R build has none).
#
# Usage:
#   Rscript plot_ld_block_audit.R <results_dir> <chr>
#
# Writes into <results_dir>:
#   share_lost_chr<chr>.png            share of local LD discarded, per
#                                      population and partition
#   boundary_max_r2_chr<chr>_<POP>.png strongest split LD at each boundary

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript plot_ld_block_audit.R <results_dir> <chr>")
out_dir <- args[[1]]
chr <- as.integer(args[[2]])

summary_df <- utils::read.csv(file.path(out_dir, sprintf("boundary_audit_summary_chr%d.csv", chr)),
  stringsAsFactors = FALSE)
boundary_df <- utils::read.csv(file.path(out_dir, sprintf("boundary_audit_per_boundary_chr%d.csv", chr)),
  stringsAsFactors = FALSE)

# Fixed categorical order (slots 1-4 of the default palette); every partition
# also gets its own point shape, so identity is never color-only.
part_names <- unique(summary_df$partition)
pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")[seq_along(part_names)]
pch <- c(16, 17, 15, 18)[seq_along(part_names)]
ink <- "#0b0b0b"
ink2 <- "#52514e"
grid_col <- "#e4e3df"
label_of <- c(fixed_250kb = "Fixed 250 kb (current)", prscs_eur = "PRS-CS EUR blocks",
  prscs_eas = "PRS-CS EAS blocks", prscs_afr = "PRS-CS AFR blocks")
labs <- ifelse(part_names %in% names(label_of), label_of[part_names], part_names)

# 1. Headline: share of local LD discarded, per population and partition.
png(file.path(out_dir, sprintf("share_lost_chr%d.png", chr)), width = 1600, height = 1000, res = 200)
par(mar = c(4.2, 5, 5, 1), bg = "#fcfcfb", fg = ink2, col.axis = ink2, col.lab = ink)
pop_order <- rev(unique(summary_df$superpop))
xmax <- max(summary_df$share_lost) * 1.1
plot(NA, xlim = c(0, xmax), ylim = c(0.5, length(pop_order) + 0.5), yaxt = "n", xaxt = "n",
  xlab = "Share of LD between SNPs within 1 Mb that is discarded", ylab = "", bty = "n")
xt <- pretty(c(0, xmax))
abline(v = xt, col = grid_col, lwd = 1)
axis(1, at = xt, labels = paste0(round(100 * xt), "%"), lwd = 0, lwd.ticks = 0)
axis(2, at = seq_along(pop_order), labels = pop_order, las = 1, lwd = 0, col.axis = ink)
offs <- seq(-0.27, 0.27, length.out = length(part_names))
for (p in seq_along(part_names)) {
  d <- summary_df[summary_df$partition == part_names[p], ]
  y <- match(d$superpop, pop_order) + offs[p]
  points(d$share_lost, y, pch = pch[p], col = pal[p], cex = 1.4)
}
title(main = sprintf("chr%d: fixed 250 kb blocks discard far more LD than published blocks", chr),
  adj = 0, line = 3.2, col.main = ink, font.main = 2, cex.main = 0.95)
legend("top", inset = c(0, -0.12), xpd = TRUE, horiz = TRUE, legend = labs, pch = pch, col = pal,
  bty = "n", text.col = ink, cex = 0.7, pt.cex = 1.1, x.intersp = 0.8, text.width = NA)
dev.off()

# 2. Per-boundary max r^2 along the chromosome, one row per partition.
xr <- range(boundary_df$bp) / 1e6
for (pop in unique(boundary_df$superpop)) {
  png(file.path(out_dir, sprintf("boundary_max_r2_chr%d_%s.png", chr, pop)), width = 1800,
    height = 300 + 330 * length(part_names), res = 200)
  par(mfrow = c(length(part_names), 1), mar = c(2.4, 4.5, 1.8, 1), oma = c(2, 0, 2, 0),
    bg = "#fcfcfb", fg = ink2, col.axis = ink2, col.lab = ink)
  for (p in seq_along(part_names)) {
    d <- boundary_df[boundary_df$superpop == pop & boundary_df$partition == part_names[p], ]
    plot(NA, xlim = xr, ylim = c(0, 1), xlab = "", ylab = "max r²", bty = "n", las = 1)
    abline(h = c(0.2, 0.5), col = grid_col)
    segments(d$bp / 1e6, 0, d$bp / 1e6, d$max_r2, col = pal[p], lwd = 1.5)
    points(d$bp / 1e6, d$max_r2, pch = 16, col = pal[p], cex = 0.6)
    title(main = sprintf("%s: %d boundaries, %d with a split pair above r² = 0.5",
      labs[p], nrow(d), sum(d$max_r2 > 0.5)), adj = 0, col.main = ink, font.main = 1, cex.main = 0.95)
  }
  mtext(sprintf("Position on chr%d (Mb, hg19)", chr), side = 1, outer = TRUE, line = 0.6, col = ink, cex = 0.8)
  mtext(sprintf("%s: strongest LD split by each boundary (SNP pairs within 1 Mb)", pop),
    side = 3, outer = TRUE, line = 0.4, adj = 0.02, col = ink, font = 2, cex = 0.9)
  dev.off()
}

message("Wrote figures to ", out_dir)
