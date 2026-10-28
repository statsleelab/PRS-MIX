#!/usr/bin/env Rscript
# Figure panel for the PRS-MIX workflow figure: a chromosome track showing
# that two jointly built PRS-CS panels share the same LD blocks, with the
# block and SNP window drawn by plot_ld_mix_heatmaps.R highlighted, so the
# track can lead into the heatmaps.
#
# Each block is drawn over the base-pair span of its SNPs (first to last
# SNP), read from each panel's HDF5 file; the script stops if the two panels'
# spans differ.
#
# Usage:
#   Rscript plot_ld_block_track.R <panel_a_h5> <panel_b_h5> <heatmap_info_tsv> <out_dir> \
#     [label_a=EUR] [label_b=ASN] [chr=22]
#
#   <heatmap_info_tsv>  ld_mix_heatmaps_info.tsv from plot_ld_mix_heatmaps.R
#                       (block, first_bp, last_bp of the window to highlight).
#
# Writes into <out_dir>: ld_block_track.pdf and ld_block_track.png.

suppressPackageStartupMessages({
  library(hdf5r)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_ld_block_track.R <panel_a_h5> <panel_b_h5> ",
       "<heatmap_info_tsv> <out_dir> [label_a=EUR] [label_b=ASN] [chr=22]", call. = FALSE)
}
arg <- function(i, default) if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
h5_a <- args[[1]]
h5_b <- args[[2]]
info <- utils::read.delim(args[[3]], stringsAsFactors = FALSE)
out_dir <- args[[4]]
label_a <- arg(5, "EUR")
label_b <- arg(6, "ASN")
chr <- arg(7, "22")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

block_spans <- function(path) {
  f <- H5File$new(path, mode = "r")
  on.exit(f$close_all())
  ids <- sort(as.integer(sub("^blk_", "", grep("^blk_[0-9]+$", names(f), value = TRUE))))
  do.call(rbind, lapply(ids, function(i) {
    bp <- f[[sprintf("blk_%d", i)]][["bp"]][]
    if (length(bp) == 0) return(NULL)
    data.frame(block = i, start = min(bp), end = max(bp), n_snps = length(bp))
  }))
}
spans_a <- block_spans(h5_a)
spans_b <- block_spans(h5_b)
if (!identical(spans_a, spans_b)) {
  stop("The two panels' block spans differ; they must be jointly built.", call. = FALSE)
}

hl_block <- info$block[[1]]
win_start <- info$first_bp[[1]]
win_end <- info$last_bp[[1]]

# ---------- Drawing ----------
# Blocks alternate two neutral tones so neighbours stay distinguishable; the
# highlighted block takes one accent colour, and the heatmap window is a dark
# bar under both tracks. Dashed lines at every block edge run through both
# tracks to show the boundaries are shared.
ink <- "#0b0b0b"
ink2 <- "#52514e"
tone <- c("#c9c7c1", "#9f9d97")
accent <- "#eda100"

rows <- c(label_a, label_b)
y_of <- c(2, 1)
names(y_of) <- rows
h <- 0.32

tracks <- do.call(rbind, lapply(rows, function(r) {
  d <- spans_a
  d$row <- r
  d$y <- y_of[[r]]
  d$fill <- ifelse(d$block == hl_block, "hl", ifelse(seq_len(nrow(d)) %% 2 == 1, "a", "b"))
  d
}))
edges <- sort(unique(c(spans_a$start, spans_a$end)))
mb <- function(x) x / 1e6

hl <- spans_a[spans_a$block == hl_block, ]
# With many blocks (e.g. ~140 250-kb tiles on chr22) the dashed edge lines
# become noise; the two identical rows already show the boundaries are
# shared, so only thin white gaps separate the tiles.
dash_edges <- nrow(spans_a) <= 40
p <- ggplot()
if (dash_edges) {
  p <- p + geom_segment(data = data.frame(x = mb(edges)),
                        aes(x = x, xend = x, y = 1 - h - 0.08, yend = 2 + h + 0.08),
                        colour = ink2, linewidth = 0.2, linetype = "22")
}
p <- p +
  geom_rect(data = tracks,
            aes(xmin = mb(start), xmax = mb(end), ymin = y - h, ymax = y + h, fill = fill),
            colour = if (dash_edges) NA else "white", linewidth = 0.1) +
  scale_fill_manual(values = c(a = tone[1], b = tone[2], hl = accent), guide = "none") +
  # Heatmap window: a dark bar under the tracks, with its label.
  annotate("rect", xmin = mb(win_start), xmax = mb(win_end),
           ymin = 1 - h - 0.30, ymax = 1 - h - 0.18, fill = ink) +
  annotate("text", x = mb((win_start + win_end) / 2), y = 1 - h - 0.42,
           label = sprintf("%d-SNP window (heatmaps)",
                           info$window_end[[1]] - info$window_start[[1]] + 1),
           size = 2.3, colour = ink, vjust = 1) +
  annotate("text", x = mb((hl$start + hl$end) / 2), y = 2 + h + 0.12,
           label = sprintf("block %d", hl_block), size = 2.3, colour = ink, vjust = 0) +
  scale_y_continuous(breaks = y_of, labels = sprintf("%s panel", names(y_of)),
                     limits = c(1 - h - 0.75, 2 + h + 0.35), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  labs(x = sprintf("Chromosome %s position (Mb)", chr), y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(colour = ink, size = 7.5, face = "bold"),
        axis.text.x = element_text(colour = ink2, size = 6.5),
        axis.title.x = element_text(colour = ink2, size = 7),
        axis.line.x = element_line(colour = ink2, linewidth = 0.3),
        axis.ticks.x = element_line(colour = ink2, linewidth = 0.3),
        plot.margin = margin(4, 6, 2, 4))

ggsave(file.path(out_dir, "ld_block_track.pdf"), p, width = 17, height = 3.4, units = "cm",
       device = grDevices::pdf)
ggsave(file.path(out_dir, "ld_block_track.png"), p, width = 17, height = 3.4, units = "cm",
       dpi = 300, bg = "white")

message(sprintf("%d shared blocks on chr%s (%.2f-%.2f Mb); highlighted block %d (%.2f-%.2f Mb)",
                nrow(spans_a), chr, mb(min(spans_a$start)), mb(max(spans_a$end)),
                hl_block, mb(hl$start), mb(hl$end)))
