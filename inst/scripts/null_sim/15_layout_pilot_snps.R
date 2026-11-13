#!/usr/bin/env Rscript
# Layout pilot, shared SNP set. Every pilot panel (each layout x panel source x
# ancestry) was built with the frozen analysis SNPs as its SNP filter. This
# step checks that each panel lists those SNPs with the analysis allele order,
# intersects them, and writes the one PRS-CS --bim_prefix every pilot arm
# uses, so the layouts differ only in their LD blocks.
#
# Usage:
#   Rscript 15_layout_pilot_snps.R <pilot_dir> <chr> <pops=EUR,AFR> <layouts=fixed250kb,fixed1000kb,eurcuts> \
#     [sources=reference,discovery]
#
# Reads <pilot_dir>/inputs/geno/analysis_snps.tsv and
# <pilot_dir>/panels/<layout>/<source>/ldblk_1kg_<pop>33kg/snpinfo_1kg_hm3; writes
# <pilot_dir>/snps/{common.bim,common.tsv,coverage.tsv}.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args("Rscript 15_layout_pilot_snps.R <pilot_dir> <chr> <pops> <layouts> [sources=reference,discovery]", 4)
pilot <- args[[1]]
chr <- as.integer(args[[2]])
pops <- toupper(strsplit(args[[3]], ",")[[1]])
layouts <- strsplit(args[[4]], ",")[[1]]
sources <- strsplit(ns_arg(args, 5, "reference,discovery"), ",")[[1]]

snps <- ns_read_tsv(file.path(pilot, "inputs", "geno", "analysis_snps.tsv"))
out <- file.path(pilot, "snps")
dir.create(out, showWarnings = FALSE)

keep <- rep(TRUE, nrow(snps))
coverage <- list()
for (layout in layouts) for (src in sources) for (pop in pops) {
  f <- file.path(ns_panel_dir(file.path(pilot, "panels", layout, src), pop), "snpinfo_1kg_hm3")
  info <- utils::read.table(f, header = TRUE, stringsAsFactors = FALSE)
  info <- info[info$CHR == chr, , drop = FALSE]
  pos <- match(info$SNP, snps$SNP)
  if (anyNA(pos)) stop(sprintf("%s lists SNPs outside the analysis set", f), call. = FALSE)
  if (any(info$A1 != snps$A1[pos] | info$A2 != snps$A2[pos] | info$BP != snps$BP[pos])) {
    stop(sprintf("%s differs from the analysis SNPs in position or allele order", f), call. = FALSE)
  }
  present <- snps$SNP %in% info$SNP
  keep <- keep & present
  coverage[[length(coverage) + 1]] <- data.frame(layout = layout, source = src, pop = pop,
                                                 n_analysis = nrow(snps), n_panel = sum(present),
                                                 n_missing = sum(!present))
}
coverage <- do.call(rbind, coverage)
common <- snps[keep, , drop = FALSE]
ns_write_tsv(coverage, file.path(out, "coverage.tsv"))
ns_write_tsv(common, file.path(out, "common.tsv"))
utils::write.table(data.frame(common$CHR, common$SNP, 0, common$BP, common$A1, common$A2),
                   file.path(out, "common.bim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print(coverage, row.names = FALSE)
message(sprintf("Common SNP set: %d of %d analysis SNPs", nrow(common), nrow(snps)))
