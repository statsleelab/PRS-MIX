#!/usr/bin/env Rscript
# Compares how independent the LD blocks of two or more built PRS-CS panels
# are: how much real LD crosses each build's block boundaries, which PRS-CS
# then treats as zero. Builds should contain the same SNPs (check with
# compare_panel_snp_sets.R) so only the block layout differs.
#
# LD comes from the HapMap3 pair files written by audit_ld_block_boundaries.R
# (SNP pairs within 1 Mb, 33KG genotypes, corrected estimator), restricted to
# the SNPs present in the built panels. Block ranges are read from each
# build's own HDF5 (start_bp/end_bp attributes), so the measurement uses
# exactly the blocks PRS-CS will see.
#
# Metrics per population and build:
#   ld_lost_pct          sum of bias-corrected r^2 over pairs split by a
#                        boundary / sum over all pairs, as a percentage
#                        (corrected r^2 = r^2 - (1 - r^2) / (n - 2)).
#   strong_split         pairs with r^2 > 0.5 placed in different blocks.
#   very_strong_split    pairs with r^2 > 0.8 placed in different blocks.
#   boundaries_over_0.5  boundaries splitting at least one pair with r^2 > 0.5.
#
# Usage:
#   Rscript compare_block_independence.R <pairs_dir> <hm3_snpinfo> <chr> <out_dir> \
#     <pops> <name1=build_dir1> <name2=build_dir2> [...]
#
#   <pairs_dir>   Directory of <POP>_w<k>.rds files from audit_ld_block_boundaries.R.
#   <hm3_snpinfo> The HapMap3 snpinfo file the pair files were indexed against.
#   <pops>        Comma-separated, e.g. EUR,ASN.
#   name=dir      Build label and the output_dir given to generate_prscs_ld_panels().
#
# Writes block_independence_summary_chr<chr>.csv, ..._boundaries_chr<chr>.csv
# (one row per boundary) and ..._blocks_chr<chr>.csv (one row per block).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript compare_block_independence.R <pairs_dir> <hm3_snpinfo> <chr> <out_dir> <pops> <name1=dir1> <name2=dir2> [...]")
}
pairs_dir <- args[[1]]
hm3_file <- args[[2]]
chr <- as.integer(args[[3]])
out_dir <- args[[4]]
pops <- toupper(strsplit(args[[5]], ",")[[1]])
builds <- strsplit(args[-(1:5)], "=", fixed = TRUE)
builds <- stats::setNames(vapply(builds, `[`, "", 2), vapply(builds, `[`, "", 1))

if (!requireNamespace("hdf5r", quietly = TRUE)) stop("Package 'hdf5r' is required.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hm3 <- utils::read.table(hm3_file, header = TRUE, stringsAsFactors = FALSE)
hm3 <- hm3[hm3$CHR == chr, , drop = FALSE]
hm3 <- hm3[order(hm3$BP), , drop = FALSE]
rownames(hm3) <- NULL
snp_key <- function(bp, a1, a2) {
  paste(bp, pmin(toupper(a1), toupper(a2)), pmax(toupper(a1), toupper(a2)))
}
hm3_key <- snp_key(hm3$BP, hm3$A1, hm3$A2)

panel_dir <- function(build, pop) file.path(build, sprintf("ldblk_1kg_%s33kg", tolower(pop)))

read_blocks <- function(build, pop) {
  h5 <- hdf5r::H5File$new(file.path(panel_dir(build, pop), sprintf("ldblk_1kg_chr%d.hdf5", chr)), mode = "r")
  on.exit(h5$close_all())
  blks <- grep("^blk_", names(h5), value = TRUE)
  blks <- blks[order(as.integer(sub("blk_", "", blks)))]
  data.frame(
    block = seq_along(blks),
    start_bp = vapply(blks, function(b) as.numeric(h5[[b]]$attr_open("start_bp")$read()), numeric(1)),
    end_bp = vapply(blks, function(b) as.numeric(h5[[b]]$attr_open("end_bp")$read()), numeric(1)),
    snps = vapply(blks, function(b) as.integer(h5[[b]]$attr_open("m_snps")$read()), integer(1)),
    row.names = NULL
  )
}

summary_rows <- list()
boundary_rows <- list()
block_rows <- list()

for (pop in pops) {
  files <- list.files(pairs_dir, pattern = sprintf("^%s_w[0-9]+\\.rds$", pop), full.names = TRUE)
  if (length(files) == 0) stop("No pair files for ", pop, " in ", pairs_dir)
  parts <- lapply(files, readRDS)
  n <- parts[[1]]$n_samples
  pairs <- do.call(rbind, lapply(parts, `[[`, "pairs"))

  # Restrict to SNPs present in the first build (builds should share SNPs).
  info <- utils::read.table(file.path(panel_dir(builds[[1]], pop), "snpinfo_1kg_hm3"),
    header = TRUE, stringsAsFactors = FALSE)
  info <- info[info$CHR == chr, , drop = FALSE]
  in_panel <- hm3_key %in% snp_key(info$BP, info$A1, info$A2)
  pairs <- pairs[in_panel[pairs$i] & in_panel[pairs$j], , drop = FALSE]
  r2 <- pairs$r^2
  r2adj <- r2 - (1 - r2) / (n - 2)
  message(sprintf("%s: n = %d, %d panel SNPs, %d pairs within 1 Mb", pop, n, sum(in_panel), nrow(pairs)))

  for (nm in names(builds)) {
    blocks <- read_blocks(builds[[nm]], pop)
    blk <- findInterval(hm3$BP, blocks$start_bp)
    bi <- blk[pairs$i]
    bj <- blk[pairs$j]
    split <- bi != bj

    # A split pair spans every boundary between its two blocks.
    big <- which(split & r2 >= 0.01)
    span <- bj[big] - bi[big]
    bidx <- rep(bi[big], span) + sequence(span) - 1L
    br2 <- rep(r2[big], span)
    nb <- nrow(blocks) - 1
    max_by <- numeric(nb)
    gt05_by <- integer(nb)
    if (length(bidx)) {
      mx <- tapply(br2, bidx, max)
      max_by[as.integer(names(mx))] <- mx
      ct <- tapply(br2 > 0.5, bidx, sum)
      gt05_by[as.integer(names(ct))] <- ct
    }

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      superpop = pop, build = nm, n_samples = n, blocks = nrow(blocks), boundaries = nb,
      snps = sum(in_panel), snps_per_block_median = stats::median(blocks$snps),
      snps_per_block_max = max(blocks$snps), pairs_within_1mb = nrow(pairs),
      ld_lost_pct = round(100 * sum(r2adj[split]) / sum(r2adj), 2),
      strong_pairs = sum(r2 > 0.5), strong_split = sum(split & r2 > 0.5),
      very_strong_pairs = sum(r2 > 0.8), very_strong_split = sum(split & r2 > 0.8),
      boundaries_over_0.5 = sum(max_by > 0.5)
    )
    boundary_rows[[length(boundary_rows) + 1]] <- data.frame(
      superpop = pop, build = nm, boundary = seq_len(nb), bp = blocks$start_bp[-1] - 0.5,
      max_r2_split = round(max_by, 4), pairs_r2_gt_0.5 = gt05_by
    )
    block_rows[[length(block_rows) + 1]] <- data.frame(superpop = pop, build = nm, blocks)
  }
}

summary_df <- do.call(rbind, summary_rows)
utils::write.csv(summary_df, file.path(out_dir, sprintf("block_independence_summary_chr%d.csv", chr)), row.names = FALSE)
utils::write.csv(do.call(rbind, boundary_rows), file.path(out_dir, sprintf("block_independence_boundaries_chr%d.csv", chr)), row.names = FALSE)
utils::write.csv(do.call(rbind, block_rows), file.path(out_dir, sprintf("block_independence_blocks_chr%d.csv", chr)), row.names = FALSE)
print(summary_df, row.names = FALSE)
message("Wrote tables to ", out_dir)
