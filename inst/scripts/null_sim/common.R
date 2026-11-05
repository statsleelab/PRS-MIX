# Shared helpers for the null-simulation scripts in this directory (see
# inst/scripts/null_simulation_plan.html). Sourced, not run.

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

# One PRS-CS output file (CHR SNP BP A1 A2 BETA, no header) as a beta vector
# over the analysis SNPs `snps`; SNPs PRS-CS dropped count as 0.
ns_read_beta <- function(f, snps) {
  d <- utils::read.table(f, header = FALSE, stringsAsFactors = FALSE,
                         col.names = c("CHR", "SNP", "BP", "A1", "A2", "BETA"))
  pos <- match(d$SNP, snps$SNP)
  if (anyNA(pos) || any(d$A1 != snps$A1[pos])) stop(sprintf("Unexpected SNPs or allele order in %s", f))
  beta <- stats::setNames(numeric(nrow(snps)), snps$SNP)
  beta[pos] <- d$BETA
  beta
}

# Posterior global shrinkage estimate from a PRS-CS log (auto runs), else NA.
ns_read_phi <- function(log_file) {
  if (!file.exists(log_file)) return(NA_real_)
  hit <- grep("global shrinkage parameter", readLines(log_file, warn = FALSE), value = TRUE)
  if (length(hit) == 0) return(NA_real_)
  as.numeric(sub(".*:\\s*([-+0-9.eE]+).*", "\\1", hit[length(hit)]))
}

# Arm name -> weight on pops[1]: mix_w<w>, single_<pops[1]> = 1, single_<pops[2]> = 0.
# A _seed<k> suffix marks a rerun of the same arm with another sampler seed.
ns_arm_weight <- function(arm, pops) {
  arm <- sub("_seed[0-9]+$", "", arm)
  if (grepl("^mix_w", arm)) as.numeric(sub("^mix_w", "", arm))
  else if (arm == paste0("single_", pops[1])) 1 else 0
}

ns_write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

ns_read_tsv <- function(path) {
  utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
}
