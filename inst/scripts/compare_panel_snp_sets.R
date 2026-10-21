#!/usr/bin/env Rscript
# Checks that two PRS-CS panel builds that differ only in block layout contain
# exactly the same variants: same SNP ids, positions, alleles, allele
# frequencies, and order. Run this before comparing PRS-CS results between the
# two builds, so any difference can be attributed to the blocks alone.
#
# Compares, for each population, the snpinfo_1kg_hm3 tables of the two builds,
# and checks each build's HDF5 blocks list the same SNPs in the same order as
# its own snpinfo table.
#
# Usage:
#   Rscript compare_panel_snp_sets.R <build_a_dir> <build_b_dir> <chr> [pops=EUR,ASN]
#
#   <build_*_dir>  output_dir passed to generate_prscs_ld_panels(), containing
#                  ldblk_1kg_<pop>33kg/ panel directories.
#
# Exits with status 1 if any check fails.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript compare_panel_snp_sets.R <build_a_dir> <build_b_dir> <chr> [pops=EUR,ASN]")
}
dir_a <- args[[1]]
dir_b <- args[[2]]
chr <- as.integer(args[[3]])
pops <- if (length(args) >= 4) toupper(strsplit(args[[4]], ",")[[1]]) else c("EUR", "ASN")

if (!requireNamespace("hdf5r", quietly = TRUE)) stop("Package 'hdf5r' is required.")

read_snpinfo <- function(build, pop) {
  d <- utils::read.table(file.path(build, sprintf("ldblk_1kg_%s33kg", tolower(pop)), "snpinfo_1kg_hm3"),
    header = TRUE, stringsAsFactors = FALSE)
  d <- d[d$CHR == chr, , drop = FALSE]
  rownames(d) <- NULL
  d
}

read_h5_snps <- function(build, pop) {
  h5 <- hdf5r::H5File$new(file.path(build, sprintf("ldblk_1kg_%s33kg", tolower(pop)),
    sprintf("ldblk_1kg_chr%d.hdf5", chr)), mode = "r")
  on.exit(h5$close_all())
  blks <- grep("^blk_", names(h5), value = TRUE)
  blks <- blks[order(as.integer(sub("blk_", "", blks)))]
  list(
    n_blocks = length(blks),
    snps = unlist(lapply(blks, function(b) h5[[b]][["snplist"]]$read())),
    bp = unlist(lapply(blks, function(b) h5[[b]][["bp"]]$read()))
  )
}

ok <- TRUE
check <- function(pass, msg) {
  cat(sprintf("[%s] %s\n", if (pass) "PASS" else "FAIL", msg))
  if (!pass) ok <<- FALSE
}

for (pop in pops) {
  a <- read_snpinfo(dir_a, pop)
  b <- read_snpinfo(dir_b, pop)
  cat(sprintf("\n%s: build A %d SNPs, build B %d SNPs\n", pop, nrow(a), nrow(b)))
  check(nrow(a) == nrow(b), sprintf("%s same number of SNPs", pop))
  for (col in c("SNP", "BP", "A1", "A2", "MAF")) {
    same <- nrow(a) == nrow(b) && identical(a[[col]], b[[col]])
    check(same, sprintf("%s identical %s column, in order", pop, col))
  }
  if (nrow(a) != nrow(b) || !identical(a$SNP, b$SNP)) {
    only_a <- setdiff(a$SNP, b$SNP)
    only_b <- setdiff(b$SNP, a$SNP)
    cat(sprintf("  only in A: %d (e.g. %s); only in B: %d (e.g. %s)\n",
      length(only_a), paste(utils::head(only_a, 3), collapse = ", "),
      length(only_b), paste(utils::head(only_b, 3), collapse = ", ")))
  }

  for (build in c(A = dir_a, B = dir_b)) {
    h <- read_h5_snps(build, pop)
    info <- if (identical(build, dir_a)) a else b
    label <- if (identical(build, dir_a)) "A" else "B"
    check(identical(h$snps, info$SNP) && identical(as.integer(h$bp), as.integer(info$BP)),
      sprintf("%s build %s: HDF5 blocks (%d) list the same SNPs, in order, as snpinfo", pop, label, h$n_blocks))
  }
}

# Panels being mixed must also match each other within a build.
for (build in c(A = dir_a, B = dir_b)) {
  label <- if (identical(build, dir_a)) "A" else "B"
  ref <- read_snpinfo(build, pops[1])
  for (pop in pops[-1]) {
    other <- read_snpinfo(build, pop)
    check(identical(ref[, c("SNP", "BP", "A1", "A2")], other[, c("SNP", "BP", "A1", "A2")]),
      sprintf("build %s: %s and %s panels hold identical SNPs (for mixing)", label, pops[1], pop))
  }
}

cat(sprintf("\nOverall: %s\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1)
