# Shared helpers for the null-simulation scripts in this directory (see
# inst/scripts/null_simulation_plan.html). Sourced, not run.

suppressPackageStartupMessages(library(gauss))

ns_args <- function(usage, n_required) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < n_required) stop(paste("Usage:", usage), call. = FALSE)
  args
}

ns_arg <- function(args, i, default) if (length(args) >= i) args[[i]] else default

ns_panel_dir <- function(panel_root, pop) file.path(panel_root, sprintf("ldblk_1kg_%s33kg", tolower(pop)))

# Reads one panel's blocks: list of (snps, ld), in block order.
ns_read_panel <- function(panel_dir, chr) {
  h5 <- hdf5r::H5File$new(file.path(panel_dir, sprintf("ldblk_1kg_chr%d.hdf5", chr)), mode = "r")
  on.exit(h5$close_all())
  blks <- grep("^blk_", names(h5), value = TRUE)
  blks <- blks[order(as.integer(sub("blk_", "", blks)))]
  lapply(blks, function(b) {
    g <- h5[[b]]
    snps <- g[["snplist"]]$read()
    ld <- if (length(snps) == 0) matrix(numeric(0), 0, 0) else g[["ldblk"]]$read()
    if (length(snps) == 1) ld <- matrix(ld, 1, 1)
    list(snps = snps, ld = ld)
  })
}

# Restricts every block to `snps` (keeping panel order) and drops empty blocks.
ns_subset_blocks <- function(blocks, snps) {
  out <- lapply(blocks, function(b) {
    keep <- which(b$snps %in% snps)
    list(snps = b$snps[keep], ld = b$ld[keep, keep, drop = FALSE])
  })
  out[vapply(out, function(b) length(b$snps) > 0, logical(1))]
}

# Dosage matrices saved by 03_extract_genotypes.R, standardized with the
# discovery set's per-ancestry mean and SD (denominator n), so the discovery
# columns satisfy sum(x^2) = n and target columns are on the same scale.
ns_load_standardized <- function(geno_dir, pop, set) {
  g <- readRDS(file.path(geno_dir, sprintf("%s_%s_dosage.rds", tolower(pop), set)))
  sc <- readRDS(file.path(geno_dir, sprintf("%s_scaling.rds", tolower(pop))))
  x <- sweep(g, 2, sc$mean, "-")
  sweep(x, 2, sc$sd, "/")
}

ns_write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

ns_read_tsv <- function(path) {
  utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
}
