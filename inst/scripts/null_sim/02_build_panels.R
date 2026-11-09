#!/usr/bin/env Rscript
# Null simulation, step 2: build aligned, HapMap3-only EUR and ASN PRS-CS
# panels from the reference set of the split only. Discovery and target
# individuals never enter the panels.
#
# Usage:
#   Rscript 02_build_panels.R <data_dir> <split_tsv> <prscs_eur_dir> <out_dir> <chr> \
#     [layout=eurcuts|fixed<N>kb] [pops=EUR,ASN] [set=reference] [snp_filter] [maf=0.01]
#
#   layout     eurcuts (PRS-CS EUR cut points) or fixed tiles of N kb, e.g.
#              fixed250kb or fixed1000kb.
#   set        Split set whose individuals build the panels. The default,
#              reference, is the only choice for analysis panels; discovery
#              builds the layout pilot's in-sample diagnostic panels.
#   snp_filter File with columns CHR BP A1 A2 (header); defaults to the
#              HapMap3 list <prscs_eur_dir>/snpinfo_1kg_hm3.
#   maf        Per-population MAF filter.
#
#   <data_dir>       Directory with 33kg_index.gz, 33kg_geno.gz, 33kg_pop_desc.txt.
#   <prscs_eur_dir>  PRS-CS's ldblk_1kg_eur, for the EUR-derived block cut points
#                    and its snpinfo_1kg_hm3 (the HapMap3 SNP filter).

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages(library(gauss))

args <- ns_args(paste(
  "Rscript 02_build_panels.R <data_dir> <split_tsv> <prscs_eur_dir> <out_dir> <chr>",
  "[layout=eurcuts|fixed<N>kb] [pops=EUR,ASN] [set=reference] [snp_filter] [maf=0.01]"
), 5)
data_dir <- args[[1]]
split <- ns_read_tsv(args[[2]])
prscs_eur <- args[[3]]
out_dir <- args[[4]]
chr <- as.integer(args[[5]])
layout <- ns_arg(args, 6, "eurcuts")
pops <- toupper(strsplit(ns_arg(args, 7, "EUR,ASN"), ",")[[1]])
set <- ns_arg(args, 8, "reference")
snp_filter <- ns_arg(args, 9, file.path(prscs_eur, "snpinfo_1kg_hm3"))
maf <- as.numeric(ns_arg(args, 10, "0.01"))
if (layout != "eurcuts" && !grepl("^fixed[0-9]+kb$", layout)) stop("layout must be eurcuts or fixed<N>kb", call. = FALSE)
block_size <- if (layout == "eurcuts") 250000 else 1000 * as.numeric(sub("^fixed([0-9]+)kb$", "\\1", layout))

keep <- sample_keep_from_split(split[split$super_population %in% pops, ], set)
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
  block_size = block_size,
  maf_cutoff = maf,
  blocks = blocks,
  snp_filter = snp_filter,
  sample_keep = keep,
  verbose = TRUE
)
message(sprintf("Built %s panels (%s, %s set: %d individuals) in %s: %s",
  layout, paste(pops, collapse = "+"), set, sum(lengths(keep)), format(Sys.time() - t0),
  paste(unlist(dirs), collapse = ", ")))
