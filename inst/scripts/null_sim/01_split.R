#!/usr/bin/env Rscript
# Null simulation, step 1: split every population of the chosen
# superpopulations into reference (50%), discovery (25%) and target (25%)
# individuals, within each population (split_reference_samples()).
#
#   reference  builds the LD panels (step 2)
#   discovery  supplies the genotypes the null GWAS is simulated on (steps 3-4)
#   target     held-out individuals scored with each arm's betas (step 6)
#
# Usage:
#   Rscript 01_split.R <pop_desc> <out_tsv> [seed=20261024] [pops=EUR,ASN]

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages(library(gauss))

args <- ns_args("Rscript 01_split.R <pop_desc> <out_tsv> [seed=20261024] [pops=EUR,ASN]", 2)
pop_desc <- args[[1]]
out_tsv <- args[[2]]
seed <- as.integer(ns_arg(args, 3, "20261024"))
pops <- toupper(strsplit(ns_arg(args, 4, "EUR,ASN"), ",")[[1]])

split <- split_reference_samples(pop_desc, pops, seed = seed)
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
ns_write_tsv(split, out_tsv)

message(sprintf("Wrote %d individuals to %s (seed %d)", nrow(split), out_tsv, seed))
print(table(split$super_population, factor(split$set, c("reference", "discovery", "target"))))
