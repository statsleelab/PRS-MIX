#!/usr/bin/env Rscript
# Validates generate_prscs_ld_panel()/generate_prscs_ld_panels()'s LD math by
# comparing it against PRS-CS's own published 1000G reference panel, for the
# same population and chromosome.
#
# Why this is the right comparison (and its limit): both panels store LD as a
# block-diagonal approximation in the same HDF5 format (blk_1..blk_N groups,
# each with an `ldblk` matrix and a `snplist`), but each panel uses its OWN
# block partition -- ours a fixed block_size tiling of the population-agnostic
# reference index, PRS-CS's published panels use variable-size, LD-structure-
# derived blocks (approximately independent LD blocks, Berisa & Pickrell 2015).
# So block INDICES do not correspond across the two files, and a naive
# blk_i-vs-blk_i comparison would be meaningless. Instead, this script pools
# every within-block SNP pair from each panel separately, then compares r for
# exactly the pairs that happen to co-occur within *some* block in BOTH panels
# -- the only pairs either file actually assigns a nonzero LD value to. The
# printed "coverage" fraction reports how much of each panel that comparison
# set covers; a low fraction is itself informative (it means the two block
# partitions disagree a lot, independent of whether the LD values that ARE
# comparable agree).
#
# Usage:
#   Rscript verify_panel_vs_prscs_reference.R \
#     <our_panel_dir> <prscs_reference_dir> <chr> <output_png> [population_label]
#
#   <our_panel_dir>       Panel directory from generate_prscs_ld_panel()/
#                         generate_prscs_ld_panels(), e.g. ldblk_1kg_eur33kg/
#                         (must contain ldblk_1kg_chr<chr>.hdf5).
#   <prscs_reference_dir> PRS-CS's own published reference panel directory for
#                         the SAME population (e.g. ldblk_1kg_eur/ from PRS-CS's
#                         1000G European reference download) -- same on-disk
#                         HDF5 format.
#   <chr>                 Chromosome to compare (one at a time; re-run per chr).
#   <output_png>           Where to write the figure.
#   [population_label]    Optional label for the figure title. Defaults to "EUR".
#
# Example:
#   Rscript inst/scripts/verify_panel_vs_prscs_reference.R \
#     prscs_ref/ldblk_1kg_eur33kg ldblk_1kg_eur 22 chr22_eur_vs_prscs.png EUR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript verify_panel_vs_prscs_reference.R <our_panel_dir> <prscs_reference_dir> <chr> <output_png> [population_label]")
}
our_panel_dir <- args[[1]]
prscs_ref_dir <- args[[2]]
chr <- as.integer(args[[3]])
output_png <- args[[4]]
population_label <- if (length(args) >= 5) args[[5]] else "EUR"

# Each panel's `snpinfo_1kg_hm3` records which allele it calls A1 for a given
# rsid -- but neither panel's own choice of A1 is guaranteed to agree with the
# other's (they're built independently, from independent per-population MAF
# computations, with no shared strand/orientation convention enforced between
# them). Reference/alt labeling is exactly the kind of thing that legitimately
# swaps between two independently-built panels for the same real variant: it
# doesn't mean the LD itself is wrong, but it DOES flip the sign of every r
# value the swapped SNP participates in relative to the other panel's
# convention. Left uncorrected, that silently cancels out real agreement when
# pooled across many SNP pairs -- so this has to be resolved before comparing
# r values at all, not treated as noise.
.read_snpinfo_alleles <- function(panel_dir, chr) {
  snpinfo <- read.table(file.path(panel_dir, "snpinfo_1kg_hm3"), header = TRUE, stringsAsFactors = FALSE)
  snpinfo <- snpinfo[snpinfo$CHR == chr, c("SNP", "A1", "A2")]
  names(snpinfo) <- c("rsid", "a1", "a2")
  snpinfo
}

if (!requireNamespace("hdf5r", quietly = TRUE)) {
  stop("Package 'hdf5r' is required. Please install it with install.packages('hdf5r').")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Please install it with install.packages('ggplot2').")
}

# Reads one HDF5 panel file into a long-format table of within-block SNP
# pairs (rsid_a < rsid_b alphabetically, so pair identity is order-independent)
# and their signed LD correlation r. Off-diagonal, cross-block pairs are not
# represented in either file's block-diagonal storage, so they simply don't
# appear here -- that's the block-partition limit described above, not a bug.
.read_panel_pairs <- function(hdf5_path) {
  h5 <- hdf5r::H5File$new(hdf5_path, mode = "r")
  on.exit(h5$close_all())

  blk_names <- grep("^blk_", names(h5), value = TRUE)
  pair_parts <- vector("list", length(blk_names))

  for (k in seq_along(blk_names)) {
    blk <- h5[[blk_names[k]]]
    snplist <- blk[["snplist"]]$read()
    m <- length(snplist)
    if (m < 2) next

    ld <- blk[["ldblk"]]$read()
    ut <- upper.tri(ld)
    idx <- which(ut, arr.ind = TRUE)

    pair_parts[[k]] <- data.frame(
      rsid_a = snplist[idx[, 1]],
      rsid_b = snplist[idx[, 2]],
      r = ld[ut],
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, pair_parts)
  if (is.null(out)) {
    return(data.frame(rsid_a = character(0), rsid_b = character(0), r = numeric(0)))
  }

  swap <- out$rsid_a > out$rsid_b
  tmp <- out$rsid_a[swap]
  out$rsid_a[swap] <- out$rsid_b[swap]
  out$rsid_b[swap] <- tmp
  out
}

message(sprintf("Reading our panel: chr%d ...", chr))
ours <- .read_panel_pairs(file.path(our_panel_dir, sprintf("ldblk_1kg_chr%d.hdf5", chr)))

message(sprintf("Reading PRS-CS reference panel: chr%d ...", chr))
theirs <- .read_panel_pairs(file.path(prscs_ref_dir, sprintf("ldblk_1kg_chr%d.hdf5", chr)))

message(sprintf("Our panel: %s within-block SNP pairs.", format(nrow(ours), big.mark = ",")))
message(sprintf("PRS-CS reference: %s within-block SNP pairs.", format(nrow(theirs), big.mark = ",")))

# Resolve each shared SNP's A1/A2 orientation between the two panels before
# comparing any r values (see rationale above .read_snpinfo_alleles).
our_alleles <- .read_snpinfo_alleles(our_panel_dir, chr)
their_alleles <- .read_snpinfo_alleles(prscs_ref_dir, chr)
allele_merge <- merge(our_alleles, their_alleles, by = "rsid", suffixes = c("_ours", "_theirs"))

concordant <- allele_merge$a1_ours == allele_merge$a1_theirs & allele_merge$a2_ours == allele_merge$a2_theirs
swapped <- allele_merge$a1_ours == allele_merge$a2_theirs & allele_merge$a2_ours == allele_merge$a1_theirs
n_unresolved <- sum(!concordant & !swapped)
if (n_unresolved > 0) {
  message(sprintf(
    "Dropping %d SNP(s) with unresolvable allele coding between panels (neither matching nor a clean swap).",
    n_unresolved
  ))
}
flip_sign <- stats::setNames(ifelse(concordant, 1, ifelse(swapped, -1, NA_real_)), allele_merge$rsid)
message(sprintf(
  "%d shared SNPs: %d concordant A1/A2, %d swapped (sign-corrected below), %d unresolved (dropped).",
  length(flip_sign), sum(concordant), sum(swapped), n_unresolved
))

merged <- merge(ours, theirs, by = c("rsid_a", "rsid_b"), suffixes = c("_ours", "_prscs"))
if (nrow(merged) == 0) {
  stop("No SNP pairs co-occur within a block in both panels -- check that both cover the same chromosome/population and use overlapping SNP sets (e.g. both restricted to HapMap3).")
}

# A pair's sign only flips if exactly one of its two SNPs is swapped between
# panels (flip_sign(a) * flip_sign(b) == -1); if both or neither are swapped,
# the pair's sign convention already agrees. Pairs touching an unresolved SNP
# (flip_sign NA) can't be sign-corrected, so they're excluded.
pair_sign <- flip_sign[merged$rsid_a] * flip_sign[merged$rsid_b]
merged$r_prscs_aligned <- merged$r_prscs * pair_sign
merged <- merged[!is.na(pair_sign), , drop = FALSE]

frac_ours <- nrow(merged) / nrow(ours)
frac_theirs <- nrow(merged) / nrow(theirs)
message(sprintf(
  "%s SNP pairs co-occur within a block in BOTH panels with resolvable allele orientation (%.1f%% of ours, %.1f%% of PRS-CS's) -- the fair comparison set.",
  format(nrow(merged), big.mark = ","), 100 * frac_ours, 100 * frac_theirs
))

r_val <- cor(merged$r_ours, merged$r_prscs_aligned)
rmse <- sqrt(mean((merged$r_ours - merged$r_prscs_aligned)^2))
message(sprintf("Pearson correlation (allele-aligned): %.4f | RMSE: %.4f", r_val, rmse))

# --- Figure --------------------------------------------------------------
# 2D binned density in a single sequential hue, not a raw scatter: at
# chromosome scale this is commonly hundreds of thousands of SNP pairs, and
# an opaque or alpha-blended scatter of that many points degrades into an
# uninformative blob long before you can see whether points track y = x.
library(ggplot2)

p <- ggplot(merged, aes(x = r_prscs_aligned, y = r_ours)) +
  geom_abline(slope = 1, intercept = 0, color = "#8a8a86", linetype = "dashed", linewidth = 0.5) +
  geom_bin2d(bins = 80) +
  scale_fill_gradientn(
    colors = c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"),
    trans = "log10",
    name = "SNP pairs"
  ) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    title = sprintf("Generated LD panel vs. PRS-CS reference panel\n(%s, chr%d)", population_label, chr),
    subtitle = sprintf(
      "r = %.3f, RMSE = %.3f, n = %s SNP pairs, allele-orientation aligned\n(%.1f%% of ours, %.1f%% of PRS-CS's)",
      r_val, rmse, format(nrow(merged), big.mark = ","), 100 * frac_ours, 100 * frac_theirs
    ),
    x = "PRS-CS reference panel r",
    y = "Our generated panel r"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "#52514e"),
    legend.position = "right"
  )

ggsave(output_png, p, width = 7.5, height = 7, dpi = 300)
message(sprintf("Wrote figure to %s", output_png))
