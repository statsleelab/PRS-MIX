#!/usr/bin/env Rscript
# Null simulation, step 2: build aligned, HapMap3-only EUR and ASN PRS-CS
# panels from the reference set of the split only. Discovery and target
# individuals never enter the panels.
#
# Usage:
#   Rscript 02_build_panels.R <data_dir> <split_tsv> <prscs_eur_dir> <out_dir> <chr> \
#     [layout=eurcuts|fixed250kb] [pops=EUR,ASN]
#
#   <data_dir>       Directory with 33kg_index.gz, 33kg_geno.gz, 33kg_pop_desc.txt.
#   <prscs_eur_dir>  PRS-CS's ldblk_1kg_eur, for the EUR-derived block cut points
#                    and its snpinfo_1kg_hm3 (the HapMap3 SNP filter).

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args(paste(
  "Rscript 02_build_panels.R <data_dir> <split_tsv> <prscs_eur_dir> <out_dir> <chr>",
  "[layout=eurcuts|fixed250kb] [pops=EUR,ASN]"
), 5)
data_dir <- args[[1]]
split <- ns_read_tsv(args[[2]])
prscs_eur <- args[[3]]
out_dir <- args[[4]]
chr <- as.integer(args[[5]])
layout <- ns_arg(args, 6, "eurcuts")
pops <- toupper(strsplit(ns_arg(args, 7, "EUR,ASN"), ",")[[1]])
if (!layout %in% c("eurcuts", "fixed250kb")) stop("layout must be eurcuts or fixed250kb", call. = FALSE)

keep <- sample_keep_from_split(split[split$super_population %in% pops, ], "reference")
blocks <- if (layout == "eurcuts") read_prscs_block_boundaries(prscs_eur, chromosomes = chr) else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
dirs <- generate_prscs_ld_panels(
  reference_index_file = file.path(data_dir, "33kg_index.gz"),
  reference_data_file = file.path(data_dir, "33kg_geno.gz"),
  reference_pop_desc_file = file.path(data_dir, "33kg_pop_desc.txt"),
  superpopulation_labels = pops,
  output_dir = out_dir,
  chromosomes = chr,
  block_size = 250000,
  maf_cutoff = 0.01,
  blocks = blocks,
  snp_filter = file.path(prscs_eur, "snpinfo_1kg_hm3"),
  sample_keep = keep,
  verbose = TRUE
)
message(sprintf("Built %s panels (%s, reference set: %d individuals) in %s: %s",
  layout, paste(pops, collapse = "+"), sum(lengths(keep)), format(Sys.time() - t0),
  paste(unlist(dirs), collapse = ", ")))
