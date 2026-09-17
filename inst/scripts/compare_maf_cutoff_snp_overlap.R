#!/usr/bin/env Rscript
# Quantifies how much maf_cutoff drives the cross-population SNP-count gap
# that generate_prscs_ld_panels() has to reconcile away.
#
# For each requested maf_cutoff, this builds each superpopulation's panel
# BOTH independently (generate_prscs_ld_panel, one call per population, no
# reconciliation) and jointly (generate_prscs_ld_panels, reconciled to a
# common SNP set per block) for the same chromosome, then reports:
#   - each population's raw SNP count at that cutoff
#   - the reconciled/common SNP count every population's aligned panel ends
#     up with
#   - what fraction of each population's raw SNPs got dropped purely by
#     reconciliation (as opposed to by the MAF filter itself)
#
# Comparing this table across cutoffs (e.g. 0.01 vs 0) shows directly whether
# lowering/removing the MAF filter closes the gap between populations, or
# whether the gap is mostly driven by SNPs that are genuinely monomorphic in
# one population regardless of cutoff.
#
# Usage:
#   Rscript compare_maf_cutoff_snp_overlap.R \
#     <33kg_index.gz> <33kg_geno.gz> <33kg_pop_desc.txt> <output_dir> \
#     <chr> <pop1,pop2,...> [maf_cutoffs] [block_size] [max_blocks_per_chr]
#
#   <maf_cutoffs>         Comma-separated list of maf_cutoff values to try.
#                         Defaults to "0.01,0".
#   <block_size>          Defaults to 250000 (package default).
#   <max_blocks_per_chr>  Optional cap, useful for a fast smoke test before
#                         committing to a full chromosome.
#
# Example:
#   Rscript inst/scripts/compare_maf_cutoff_snp_overlap.R \
#     /shared/leed13_shared/1KG_LD/33kg_index.gz \
#     /shared/leed13_shared/1KG_LD/33kg_geno.gz \
#     /shared/leed13_shared/1KG_LD/33kg_pop_desc.txt \
#     /tmp/maf_cutoff_compare 22 EUR,ASN

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript compare_maf_cutoff_snp_overlap.R <33kg_index.gz> <33kg_geno.gz> <33kg_pop_desc.txt> <output_dir> <chr> <pop1,pop2,...> [maf_cutoffs] [block_size] [max_blocks_per_chr]")
}

index_file <- args[[1]]
geno_file <- args[[2]]
pop_desc_file <- args[[3]]
output_dir <- args[[4]]
chr <- as.integer(args[[5]])
pops <- strsplit(args[[6]], ",")[[1]]
maf_cutoffs <- if (length(args) >= 7) as.numeric(strsplit(args[[7]], ",")[[1]]) else c(0.01, 0)
block_size <- if (length(args) >= 8) as.numeric(args[[8]]) else 250000
max_blocks_per_chr <- if (length(args) >= 9) as.integer(args[[9]]) else NULL

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_snpinfo_count <- function(panel_dir) {
  path <- file.path(panel_dir, "snpinfo_1kg_hm3")
  nrow(utils::read.table(path, header = TRUE, stringsAsFactors = FALSE))
}

results <- list()

for (cutoff in maf_cutoffs) {
  message(sprintf("\n===== maf_cutoff = %s =====", cutoff))
  cutoff_dir <- file.path(output_dir, sprintf("maf_%s", gsub("\\.", "p", as.character(cutoff))))

  raw_dir <- file.path(cutoff_dir, "raw")
  raw_n <- stats::setNames(numeric(length(pops)), pops)
  for (pop in pops) {
    message(sprintf("-- raw (unaligned) panel: %s --", pop))
    panel_dir <- gauss::generate_prscs_ld_panel(
      reference_index_file = index_file,
      reference_data_file = geno_file,
      reference_pop_desc_file = pop_desc_file,
      superpopulation_label = pop,
      output_dir = raw_dir,
      chromosomes = chr,
      block_size = block_size,
      maf_cutoff = cutoff,
      max_blocks_per_chr = max_blocks_per_chr,
      verbose = TRUE
    )
    raw_n[[pop]] <- read_snpinfo_count(panel_dir)
  }

  message("-- aligned (reconciled) panels --")
  aligned_dir <- file.path(cutoff_dir, "aligned")
  panel_dirs <- gauss::generate_prscs_ld_panels(
    reference_index_file = index_file,
    reference_data_file = geno_file,
    reference_pop_desc_file = pop_desc_file,
    superpopulation_labels = pops,
    output_dir = aligned_dir,
    chromosomes = chr,
    block_size = block_size,
    maf_cutoff = cutoff,
    max_blocks_per_chr = max_blocks_per_chr,
    verbose = TRUE
  )
  # Reconciliation guarantees every population's aligned panel has the same
  # SNP count, so reading the first is representative of all of them.
  common_n <- read_snpinfo_count(panel_dirs[[1]])

  results[[as.character(cutoff)]] <- data.frame(
    maf_cutoff = cutoff,
    population = pops,
    raw_snps = as.integer(raw_n[pops]),
    common_snps = as.integer(common_n),
    pct_dropped_by_reconciliation = round(100 * (1 - common_n / raw_n[pops]), 2),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, results)
rownames(summary_table) <- NULL

message("\n===== Summary: chr", chr, ", populations ", paste(pops, collapse = ", "), " =====")
print(summary_table)

out_csv <- file.path(output_dir, "maf_cutoff_snp_overlap_summary.csv")
utils::write.csv(summary_table, out_csv, row.names = FALSE)
message("\nWrote summary table to ", out_csv)
