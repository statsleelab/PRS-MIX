#!/usr/bin/env Rscript
# Summarizes the LD block-boundary audit: given the per-window HapMap3 pair LD
# written by audit_ld_block_boundaries.R, measures how much LD each candidate
# block partition discards at its boundaries, for each superpopulation.
#
# Partitions compared (each applied identically to every population, since
# parse_ldblk_mix needs shared blocks):
#   fixed_250kb  -- the current generate_prscs_ld_panels() tiling: 250 kb tiles
#                   starting at the chromosome's first bp in the 33KG index.
#   prscs_<pop>  -- the LD-derived blocks in PRS-CS's published 1000 Genomes
#                   panels (ldblk_1kg_eur / _eas / _afr), recovered from each
#                   block's SNP positions; a boundary is the midpoint between
#                   the last SNP of one block and the first SNP of the next.
#
# Metrics, all over HapMap3 SNP pairs at most step_bp apart:
#   share_lost        sum of bias-corrected r^2 over pairs split by a boundary
#                     / sum over all pairs. Corrected r^2 = r^2 - (1-r^2)/(n-2)
#                     removes the ~1/n upward bias of sample r^2, which would
#                     otherwise inflate sums in small populations.
#   snp_frac_*        per-SNP version: that SNP's share of local LD lying
#                     outside its block (SNPs with local corrected sum >= 0.2).
#   n_pairs_gt_*      split pairs with raw r^2 above 0.1 / 0.2 / 0.5.
#   max_r2_split      largest raw r^2 among split pairs.
#   blocks, SNPs per block, and relative storage (sum m^2) / sampler cost
#   (sum m^3) for the computational trade-off.
#
# Writes CSV tables only; plot_ld_block_audit.R draws the figures (it can run
# anywhere, including machines whose R lacks a PNG device, as RedHawk's does).
#
# Usage:
#   Rscript summarize_ld_block_audit.R <pairs_dir> <hm3_snpinfo> <chr> \
#     <prscs_ref_root> <index_chr_min_bp> <out_dir> [pops=EUR,ASN,AFR,SAS] \
#     [block_size=250000]
#
#   <pairs_dir>         Directory of <POP>_w<k>.rds files from the worker.
#   <prscs_ref_root>    Directory containing ldblk_1kg_eur/, ldblk_1kg_eas/,
#                       ldblk_1kg_afr/ (PRS-CS's published panels).
#   <index_chr_min_bp>  First bp of this chromosome in the 33KG index (where
#                       the fixed tiling starts).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript summarize_ld_block_audit.R <pairs_dir> <hm3_snpinfo> <chr> <prscs_ref_root> <index_chr_min_bp> <out_dir> [pops=EUR,ASN,AFR,SAS] [block_size=250000]")
}
pairs_dir <- args[[1]]
hm3_file <- args[[2]]
chr <- as.integer(args[[3]])
prscs_root <- args[[4]]
index_min_bp <- as.numeric(args[[5]])
out_dir <- args[[6]]
pops <- if (length(args) >= 7) toupper(strsplit(args[[7]], ",")[[1]]) else c("EUR", "ASN", "AFR", "SAS")
block_size <- if (length(args) >= 8) as.numeric(args[[8]]) else 250000

if (!requireNamespace("hdf5r", quietly = TRUE)) stop("Package 'hdf5r' is required.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hm3 <- utils::read.table(hm3_file, header = TRUE, stringsAsFactors = FALSE)
hm3 <- hm3[hm3$CHR == chr, , drop = FALSE]
hm3 <- hm3[order(hm3$BP), , drop = FALSE]
rownames(hm3) <- NULL

# ---- Partitions: each is a sorted vector of boundary positions ------------

prscs_boundaries <- function(ref_dir) {
  info <- utils::read.table(file.path(ref_dir, "snpinfo_1kg_hm3"), header = TRUE,
    stringsAsFactors = FALSE)
  info <- info[info$CHR == chr, , drop = FALSE]
  bp_by_snp <- stats::setNames(info$BP, info$SNP)
  h5 <- hdf5r::H5File$new(file.path(ref_dir, sprintf("ldblk_1kg_chr%d.hdf5", chr)), mode = "r")
  on.exit(h5$close_all())
  blks <- grep("^blk_", names(h5), value = TRUE)
  blks <- blks[order(as.integer(sub("blk_", "", blks)))]
  rng <- t(vapply(blks, function(b) {
    bp <- bp_by_snp[h5[[b]][["snplist"]]$read()]
    bp <- bp[!is.na(bp)]
    if (length(bp) == 0) c(NA_real_, NA_real_) else range(bp)
  }, numeric(2)))
  rng <- rng[!is.na(rng[, 1]), , drop = FALSE]
  if (is.unsorted(rng[, 1]) || any(rng[-1, 1] <= rng[-nrow(rng), 2])) {
    stop(sprintf("Blocks in %s are not ordered and non-overlapping.", ref_dir))
  }
  (rng[-nrow(rng), 2] + rng[-1, 1]) / 2
}

fixed_boundaries <- function() {
  b <- seq(index_min_bp, max(hm3$BP) + block_size, by = block_size) - 0.5
  b[b > min(hm3$BP) & b < max(hm3$BP)]
}

partitions <- list(fixed_250kb = fixed_boundaries())
for (p in c("eur", "eas", "afr")) {
  ref_dir <- file.path(prscs_root, sprintf("ldblk_1kg_%s", p))
  if (dir.exists(ref_dir)) partitions[[paste0("prscs_", p)]] <- prscs_boundaries(ref_dir)
}
if (block_size != 250000) names(partitions)[1] <- sprintf("fixed_%gkb", block_size / 1000)
message("Partitions: ", paste(sprintf("%s (%d boundaries)", names(partitions),
  vapply(partitions, length, integer(1))), collapse = ", "))

# Boundary list, for reference and for the per-boundary table.
utils::write.csv(do.call(rbind, lapply(names(partitions), function(nm) {
  data.frame(partition = nm, boundary = seq_along(partitions[[nm]]), bp = partitions[[nm]])
})), file.path(out_dir, sprintf("partitions_chr%d.csv", chr)), row.names = FALSE)

# ---- Per-population metrics ----------------------------------------------

# Sum of w by HapMap3 row index, as a full-length vector.
snp_sum <- function(idx, w) {
  out <- numeric(nrow(hm3))
  if (length(idx)) {
    s <- rowsum(w, idx)
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}

summary_rows <- list()
boundary_rows <- list()

for (pop in pops) {
  files <- list.files(pairs_dir, pattern = sprintf("^%s_w[0-9]+\\.rds$", pop), full.names = TRUE)
  if (length(files) == 0) {
    warning("No pair files for ", pop, "; skipping.")
    next
  }
  parts <- lapply(files, readRDS)
  n <- parts[[1]]$n_samples
  pairs <- do.call(rbind, lapply(parts, `[[`, "pairs"))
  present <- sort(unique(unlist(lapply(parts, `[[`, "present"))))
  message(sprintf("%s: n = %d, %d HapMap3 SNPs present, %d pairs", pop, n, length(present), nrow(pairs)))

  r2 <- pairs$r^2
  r2adj <- r2 - (1 - r2) / (n - 2)
  local_total <- snp_sum(pairs$i, r2adj) + snp_sum(pairs$j, r2adj)

  for (nm in names(partitions)) {
    bnd <- partitions[[nm]]
    blk <- findInterval(hm3$BP, bnd) + 1L
    bi <- blk[pairs$i]
    bj <- blk[pairs$j]
    split_pair <- bi != bj

    out_snp <- snp_sum(pairs$i[split_pair], r2adj[split_pair]) +
      snp_sum(pairs$j[split_pair], r2adj[split_pair])
    eligible <- present[local_total[present] >= 0.2]
    frac <- out_snp[eligible] / local_total[eligible]

    m <- tabulate(blk[present])
    m <- m[m > 0]

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      superpop = pop, n_samples = n, partition = nm,
      n_boundaries = length(bnd), n_blocks = length(m),
      snps_per_block_median = stats::median(m), snps_per_block_max = max(m),
      sum_m2 = sum(as.numeric(m)^2), sum_m3 = sum(as.numeric(m)^3),
      share_lost = sum(r2adj[split_pair]) / sum(r2adj),
      snp_frac_median = stats::median(frac),
      snp_frac_p90 = unname(stats::quantile(frac, 0.9)),
      snp_n_eligible = length(eligible),
      n_pairs_gt_0.1 = sum(r2[split_pair] > 0.1),
      n_pairs_gt_0.2 = sum(r2[split_pair] > 0.2),
      n_pairs_gt_0.5 = sum(r2[split_pair] > 0.5),
      max_r2_split = if (any(split_pair)) max(r2[split_pair]) else 0
    )

    # Per boundary: a split pair spans every boundary between its two blocks.
    # Only pairs with r^2 >= 0.01 are expanded; boundaries with none report 0.
    big <- which(split_pair & r2 >= 0.01)
    span <- bj[big] - bi[big]
    bidx <- rep(bi[big], span) + sequence(span) - 1L
    br2 <- rep(r2[big], span)
    max_by <- numeric(length(bnd))
    gt02_by <- integer(length(bnd))
    if (length(bidx)) {
      mx <- tapply(br2, bidx, max)
      max_by[as.integer(names(mx))] <- mx
      ct <- tapply(br2 > 0.2, bidx, sum)
      gt02_by[as.integer(names(ct))] <- ct
    }
    boundary_rows[[length(boundary_rows) + 1]] <- data.frame(
      superpop = pop, partition = nm, boundary = seq_along(bnd), bp = bnd,
      max_r2 = max_by, n_pairs_gt_0.2 = gt02_by
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
boundary_df <- do.call(rbind, boundary_rows)
base_cost <- summary_df[summary_df$partition == names(partitions)[1], c("superpop", "sum_m2", "sum_m3")]
summary_df$storage_vs_fixed <- summary_df$sum_m2 / base_cost$sum_m2[match(summary_df$superpop, base_cost$superpop)]
summary_df$cost_vs_fixed <- summary_df$sum_m3 / base_cost$sum_m3[match(summary_df$superpop, base_cost$superpop)]
summary_df$sum_m2 <- NULL
summary_df$sum_m3 <- NULL

utils::write.csv(summary_df, file.path(out_dir, sprintf("boundary_audit_summary_chr%d.csv", chr)), row.names = FALSE)
utils::write.csv(boundary_df, file.path(out_dir, sprintf("boundary_audit_per_boundary_chr%d.csv", chr)), row.names = FALSE)
worst <- do.call(rbind, lapply(split(boundary_df, list(boundary_df$superpop, boundary_df$partition), drop = TRUE),
  function(d) utils::head(d[order(-d$max_r2), ], 10)))
utils::write.csv(worst, file.path(out_dir, sprintf("boundary_audit_worst_chr%d.csv", chr)), row.names = FALSE)
print(summary_df, row.names = FALSE, digits = 3)

message("Wrote tables to ", out_dir, "; make figures with plot_ld_block_audit.R")
