#!/usr/bin/env Rscript
# Figure panels for the PRS-MIX workflow figure: signed LD heatmaps for one
# window of one block in two aligned ancestry panels, and their weighted mix,
#   w * D_EUR + (1 - w) * D_ASN = D_mix
# drawn from PRS-CS-format panels written by generate_prscs_ld_panels(). The
# panels must be jointly built, so every block has the same SNPs in the same
# order in both (checked below).
#
# Usage:
#   Rscript plot_ld_mix_heatmaps.R <panel_a_h5> <panel_b_h5> <out_dir> \
#     [w_a=0.7] [label_a=EUR] [label_b=ASN] [block=auto] [window=150]
#
#   <panel_a_h5>, <panel_b_h5>  ldblk_1kg_chr<chr>.hdf5 from each panel directory.
#   w_a      mixing proportion of panel A (panel B gets 1 - w_a).
#   block    block index, or "auto" for the block whose best window differs
#            most between the two panels (mean |r_A - r_B|).
#   window   number of consecutive SNPs to show; the window within the block
#            with the largest mean |r_A - r_B| is used. Use 0 for the whole block.
#
# Writes into <out_dir>: ld_mix_heatmaps.pdf and ld_mix_heatmaps.png, plus
# ld_mix_heatmaps_info.tsv recording the block, window and SNP range drawn.

suppressPackageStartupMessages({
  library(hdf5r)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_ld_mix_heatmaps.R <panel_a_h5> <panel_b_h5> <out_dir> ",
       "[w_a=0.7] [label_a=EUR] [label_b=ASN] [block=auto] [window=150]", call. = FALSE)
}
arg <- function(i, default) if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
h5_a <- args[[1]]
h5_b <- args[[2]]
out_dir <- args[[3]]
w_a <- as.numeric(arg(4, "0.7"))
label_a <- arg(5, "EUR")
label_b <- arg(6, "ASN")
block_arg <- arg(7, "auto")
window <- as.integer(arg(8, "150"))
if (!is.finite(w_a) || w_a < 0 || w_a > 1) stop("w_a must be in [0, 1]", call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fa <- H5File$new(h5_a, mode = "r")
fb <- H5File$new(h5_b, mode = "r")
on.exit({ fa$close_all(); fb$close_all() }, add = TRUE)

n_blk <- length(grep("^blk_[0-9]+$", names(fa)))
if (n_blk != length(grep("^blk_[0-9]+$", names(fb)))) {
  stop("The two panels have different numbers of blocks; they must be jointly built.", call. = FALSE)
}

read_block <- function(i) {
  g <- sprintf("blk_%d", i)
  snps_a <- fa[[g]][["snplist"]][]
  snps_b <- fb[[g]][["snplist"]][]
  if (!identical(snps_a, snps_b)) {
    stop(sprintf("Block %d has different SNP lists in the two panels.", i), call. = FALSE)
  }
  if (length(snps_a) < 2) return(NULL)
  list(snps = snps_a, bp = fa[[g]][["bp"]][],
       ra = fa[[g]][["ldblk"]][, ], rb = fb[[g]][["ldblk"]][, ])
}

# Best window (start index and score) within a block: largest mean |r_A - r_B|
# over pairs inside the window.
best_window <- function(blk) {
  m <- length(blk$snps)
  if (window <= 0 || window >= m) {
    ut <- upper.tri(blk$ra)
    return(c(start = 1, end = m, score = mean(abs(blk$ra[ut] - blk$rb[ut]))))
  }
  d <- abs(blk$ra - blk$rb)
  starts <- seq(1, m - window + 1, by = max(1, window %/% 10))
  scores <- vapply(starts, function(s) {
    idx <- s:(s + window - 1)
    sub <- d[idx, idx]
    mean(sub[upper.tri(sub)])
  }, numeric(1))
  k <- which.max(scores)
  c(start = starts[k], end = starts[k] + window - 1, score = scores[k])
}

if (identical(block_arg, "auto")) {
  cand <- lapply(seq_len(n_blk), function(i) {
    blk <- read_block(i)
    if (is.null(blk)) return(NULL)
    c(block = i, best_window(blk))
  })
  cand <- do.call(rbind, cand)
  block <- as.integer(cand[which.max(cand[, "score"]), "block"])
} else {
  block <- as.integer(block_arg)
}
blk <- read_block(block)
if (is.null(blk)) stop(sprintf("Block %d has fewer than two SNPs.", block), call. = FALSE)
win <- best_window(blk)
idx <- win[["start"]]:win[["end"]]

ra <- blk$ra[idx, idx]
rb <- blk$rb[idx, idx]
rmix <- w_a * ra + (1 - w_a) * rb

info <- data.frame(
  block = block, n_snps_block = length(blk$snps),
  window_start = win[["start"]], window_end = win[["end"]],
  first_snp = blk$snps[idx[1]], last_snp = blk$snps[idx[length(idx)]],
  first_bp = blk$bp[idx[1]], last_bp = blk$bp[idx[length(idx)]],
  mean_abs_diff = win[["score"]], w_a = w_a
)
utils::write.table(info, file.path(out_dir, "ld_mix_heatmaps_info.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- Drawing ----------
# Signed r on a diverging scale: blue (negative) - neutral gray (0) - red
# (positive), equal arms, fixed limits so all three panels share one legend.
ink <- "#0b0b0b"
ink2 <- "#52514e"
div_low <- "#104281"
div_mid <- "#f0efec"
div_high <- "#a3201b"

# Lower triangle only (including the diagonal): the matrices are symmetric,
# so the upper triangle would only repeat it.
to_long <- function(r) {
  n <- nrow(r)
  d <- data.frame(i = rep(seq_len(n), times = n), j = rep(seq_len(n), each = n),
                  r = as.vector(r))
  d[d$j >= d$i, ]
}

heat <- function(r, title) {
  ggplot(to_long(r), aes(i, j, fill = r)) +
    geom_raster() +
    scale_fill_gradient2(low = div_low, mid = div_mid, high = div_high, midpoint = 0,
                         limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                         name = expression(LD ~ italic(r))) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_fixed() +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(hjust = 0.5, colour = ink, face = "bold", size = 9,
                                    margin = margin(b = 3)),
          legend.title = element_text(colour = ink2, size = 8),
          legend.text = element_text(colour = ink2, size = 7))
}

op <- function(sym) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = sym, size = 3.6, colour = ink) +
    coord_cartesian(clip = "off") +
    theme_void()
}

fmt_w <- function(w) formatC(w, format = "fg", digits = 2)
# The weights sit between the panels ("0.7 x EUR + 0.3 x ASN = mixed"); each
# heatmap shows its panel's unscaled LD.
p <- op(sprintf("%s ×", fmt_w(w_a))) + heat(ra, label_a) +
  op(sprintf("+  %s ×", fmt_w(1 - w_a))) + heat(rb, label_b) + op("=") +
  heat(rmix, "Mixed LD") +
  plot_layout(nrow = 1, widths = c(0.3, 1, 0.45, 1, 0.15, 1),
              guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "ld_mix_heatmaps.pdf"), p, width = 17, height = 5.2, units = "cm",
       device = grDevices::pdf)
ggsave(file.path(out_dir, "ld_mix_heatmaps.png"), p, width = 17, height = 5.2, units = "cm",
       dpi = 300, bg = "white")

message(sprintf("Block %d (%d SNPs), SNPs %d-%d (bp %d-%d), mean |r_%s - r_%s| = %.3f",
                block, length(blk$snps), win[["start"]], win[["end"]],
                blk$bp[idx[1]], blk$bp[idx[length(idx)]], label_a, label_b, win[["score"]]))
