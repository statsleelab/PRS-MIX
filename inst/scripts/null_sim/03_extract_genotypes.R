#!/usr/bin/env Rscript
# Null simulation, step 3: read the discovery and target individuals' dosages
# at the panels' SNPs, fix the analysis SNP set, and record each ancestry's
# discovery scaling constants.
#
# The analysis SNP set is the panels' (identical, aligned) SNP set minus SNPs
# with any missing call in a discovery or target set, or no variation in a
# discovery set. Every later step -- simulated sumstats, every PRS-CS arm,
# LD-fit scores and target scores -- uses exactly this set.
#
# Standardization (plan section 5.1): per ancestry, center each SNP at the
# discovery set's mean dosage and divide by its SD with denominator n, so
# discovery columns satisfy sum(x^2) = n. Target individuals use the same
# constants.
#
# Also writes a preflight comparison of the panel estimator (CalWgtCov on the
# discovery set) against the discovery Gram correlation Z'Z/n for a few blocks
# (plan section 4.3): these differ slightly in finite samples by design.
#
# Usage:
#   Rscript 03_extract_genotypes.R <data_dir> <split_tsv> <panel_root> <chr> <out_dir> [pops=EUR,ASN]
#
# Outputs in <out_dir>: <pop>_{discovery,target}_dosage.rds (integer matrices,
# individuals x analysis SNPs), <pop>_{discovery,target}_samples.tsv,
# <pop>_scaling.rds, analysis_snps.tsv, analysis_snps.bim (PRS-CS
# --bim_prefix), excluded_snps.tsv, preflight_ld_estimator.tsv.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args("Rscript 03_extract_genotypes.R <data_dir> <split_tsv> <panel_root> <chr> <out_dir> [pops=EUR,ASN]", 5)
data_dir <- args[[1]]
split <- ns_read_tsv(args[[2]])
panel_root <- args[[3]]
chr <- as.integer(args[[4]])
out_dir <- args[[5]]
pops <- toupper(strsplit(ns_arg(args, 6, "EUR,ASN"), ",")[[1]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ref_files <- list(
  reference_index_file = file.path(data_dir, "33kg_index.gz"),
  reference_data_file = file.path(data_dir, "33kg_geno.gz"),
  reference_pop_desc_file = file.path(data_dir, "33kg_pop_desc.txt")
)

read_info <- function(pop) {
  info <- utils::read.table(file.path(ns_panel_dir(panel_root, pop), "snpinfo_1kg_hm3"),
                            header = TRUE, stringsAsFactors = FALSE)
  info[info$CHR == chr, , drop = FALSE]
}
info <- read_info(pops[1])
for (pop in pops[-1]) {
  other <- read_info(pop)
  if (!identical(other[, c("SNP", "BP", "A1", "A2")], info[, c("SNP", "BP", "A1", "A2")])) {
    stop(sprintf("Panels for %s and %s do not list identical SNPs; rebuild them together.", pops[1], pop),
         call. = FALSE)
  }
}
info$key <- paste(info$BP, info$A1, info$A2, sep = ":")
message(sprintf("chr%d: %d panel SNPs", chr, nrow(info)))

read_set <- function(pop, set) {
  keep <- sample_keep_from_split(split[split$super_population == pop, ], set)
  g <- do.call(readSuperPopGenotypes, c(list(
    chr = chr, start_bp = min(info$BP), end_bp = max(info$BP),
    superpopulation_label = pop, keep_bp = info$BP, sample_keep = keep
  ), ref_files))
  key <- paste(g$snplist$bp, g$snplist$a1, g$snplist$a2, sep = ":")
  pos <- match(info$key, key)
  if (anyNA(pos)) stop(sprintf("%d panel SNPs not found in the reference index for %s.", sum(is.na(pos)), pop))
  list(dosage = g$genotypes[, pos, drop = FALSE], samples = g$samples)
}

sets <- list()
excl <- list()
for (pop in pops) {
  for (set in c("discovery", "target")) {
    message(sprintf("Reading %s %s genotypes...", pop, set))
    s <- read_set(pop, set)
    sets[[pop]][[set]] <- s
    has_na <- colSums(is.na(s$dosage)) > 0
    if (any(has_na)) excl[[length(excl) + 1]] <- data.frame(SNP = info$SNP[has_na], reason = sprintf("missing_%s_%s", pop, set))
  }
  sd0 <- apply(sets[[pop]]$discovery$dosage, 2, function(x) stats::sd(x, na.rm = TRUE))
  mono <- !is.na(sd0) & sd0 == 0
  if (any(mono)) excl[[length(excl) + 1]] <- data.frame(SNP = info$SNP[mono], reason = sprintf("monomorphic_%s_discovery", pop))
}
excl <- if (length(excl)) do.call(rbind, excl) else data.frame(SNP = character(0), reason = character(0))
keep_snp <- !info$SNP %in% excl$SNP
ns_write_tsv(excl, file.path(out_dir, "excluded_snps.tsv"))
message(sprintf("Analysis SNP set: %d of %d panel SNPs (%d excluded)", sum(keep_snp), nrow(info), sum(!keep_snp)))

analysis <- info[keep_snp, c("CHR", "SNP", "BP", "A1", "A2"), drop = FALSE]
ns_write_tsv(analysis, file.path(out_dir, "analysis_snps.tsv"))
utils::write.table(data.frame(analysis$CHR, analysis$SNP, 0, analysis$BP, analysis$A1, analysis$A2),
                   file.path(out_dir, "analysis_snps.bim"), sep = "\t", quote = FALSE,
                   row.names = FALSE, col.names = FALSE)

for (pop in pops) {
  disc <- sets[[pop]]$discovery$dosage[, keep_snp, drop = FALSE]
  n <- nrow(disc)
  mu <- colMeans(disc)
  sdn <- sqrt(colMeans(sweep(disc, 2, mu)^2))
  saveRDS(list(mean = mu, sd = sdn, n = n, snps = analysis$SNP), file.path(out_dir, sprintf("%s_scaling.rds", tolower(pop))))
  for (set in c("discovery", "target")) {
    d <- sets[[pop]][[set]]$dosage[, keep_snp, drop = FALSE]
    colnames(d) <- analysis$SNP
    saveRDS(d, file.path(out_dir, sprintf("%s_%s_dosage.rds", tolower(pop), set)), compress = FALSE)
    ns_write_tsv(sets[[pop]][[set]]$samples, file.path(out_dir, sprintf("%s_%s_samples.tsv", tolower(pop), set)))
  }
  message(sprintf("%s: %d discovery, %d target individuals", pop, n, nrow(sets[[pop]]$target$dosage)))
}

# Preflight: panel estimator vs. Gram correlation on the same discovery set,
# for the first, middle and last non-trivial panel blocks.
blocks <- ns_subset_blocks(ns_read_panel(ns_panel_dir(panel_root, pops[1]), chr), analysis$SNP)
sizes <- vapply(blocks, function(b) length(b$snps), integer(1))
pick <- unique(which(sizes >= 10)[c(1, ceiling(sum(sizes >= 10) / 2), sum(sizes >= 10))])
pre <- list()
for (pop in pops) {
  Z <- ns_load_standardized(out_dir, pop, "discovery")
  keep <- sample_keep_from_split(split[split$super_population == pop, ], "discovery")
  for (b in pick) {
    snps <- blocks[[b]]$snps
    bp <- analysis$BP[match(snps, analysis$SNP)]
    est <- do.call(computeSuperPopLDWindow, c(list(
      chr = chr, start_bp = min(bp), end_bp = max(bp), superpopulation_label = pop,
      maf_cutoff = 0, keep_bp = bp, min_snps = 1L, sample_keep = keep
    ), ref_files))
    pos <- match(snps, est$snplist$rsid)
    ok <- !is.na(pos)
    gram <- crossprod(Z[, snps[ok], drop = FALSE]) / nrow(Z)
    d <- abs(est$cormat[pos[ok], pos[ok]] - gram)
    pre[[length(pre) + 1]] <- data.frame(pop = pop, block = b, m = sum(ok),
                                         max_abs_diff = max(d), mean_abs_diff = mean(d[upper.tri(d)]))
  }
}
pre <- do.call(rbind, pre)
ns_write_tsv(pre, file.path(out_dir, "preflight_ld_estimator.tsv"))
print(pre)
