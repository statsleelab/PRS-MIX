#!/usr/bin/env Rscript
# Validates a set of aligned PRS-CS panels built together by
# generate_prscs_ld_panels() for one chromosome, one population at a time:
#
#   1. Structure: every population's panel lists the same SNPs, alleles and
#      order, and each panel's HDF5 blocks match its own snpinfo file.
#   2. Matrices: every block's LD matrix is finite, symmetric, has a unit
#      diagonal, has off-diagonal |r| <= 1, and is positive definite (PRS-CS
#      takes a Cholesky factor of each block).
#   3. Agreement with PRS-CS's published 1000 Genomes panel for the matching
#      population, where one exists: correlation and RMSE of r over SNP pairs
#      that share a block in both panels, after aligning allele orientation.
#   4. Population contrasts: LD between neighbouring SNPs and the allele
#      frequency spectrum, which should differ between populations in the
#      expected direction (e.g. weaker short-range LD in AFR).
#
# Usage:
#   Rscript validate_prscs_panels.R <build_dir> <chr> <pops> [prscs_root] [ref_map] [pop_desc]
#
#   <build_dir>   output_dir given to generate_prscs_ld_panels().
#   <pops>        Comma-separated superpopulations, e.g. EUR,ASN,AFR,AMR,SAS.
#   [prscs_root]  Directory holding PRS-CS's ldblk_1kg_<pop>/ panels.
#   [ref_map]     Our population -> PRS-CS panel, defaults to EUR=eur,ASN=eas,AFR=afr.
#   [pop_desc]    33KG population description file, to report sample sizes.
#
# Exits with status 1 if a structural or matrix check fails.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript validate_prscs_panels.R <build_dir> <chr> <pops> [prscs_root] [ref_map] [pop_desc]")
}
build_dir <- args[[1]]
chr <- as.integer(args[[2]])
pops <- toupper(strsplit(args[[3]], ",")[[1]])
prscs_root <- if (length(args) >= 4) args[[4]] else NA
ref_map <- if (length(args) >= 5) args[[5]] else "EUR=eur,ASN=eas,AFR=afr"
ref_map <- strsplit(strsplit(ref_map, ",")[[1]], "=")
ref_map <- stats::setNames(vapply(ref_map, `[`, "", 2), toupper(vapply(ref_map, `[`, "", 1)))
pop_desc_file <- if (length(args) >= 6) args[[6]] else NA

if (!requireNamespace("hdf5r", quietly = TRUE)) stop("Package 'hdf5r' is required.")

ok <- TRUE
check <- function(pass, msg) {
  cat(sprintf("  [%s] %s\n", if (pass) "PASS" else "FAIL", msg))
  if (!isTRUE(pass)) ok <<- FALSE
}

read_panel <- function(dir, tag = "1kg") {
  h5 <- hdf5r::H5File$new(file.path(dir, sprintf("ldblk_%s_chr%d.hdf5", tag, chr)), mode = "r")
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

read_snpinfo <- function(dir, tag = "1kg") {
  d <- utils::read.table(file.path(dir, sprintf("snpinfo_%s_hm3", tag)), header = TRUE, stringsAsFactors = FALSE)
  d[d$CHR == chr, , drop = FALSE]
}

panel_dir <- function(pop) file.path(build_dir, sprintf("ldblk_1kg_%s33kg", tolower(pop)))

panels <- lapply(stats::setNames(pops, pops), function(p) read_panel(panel_dir(p)))
infos <- lapply(stats::setNames(pops, pops), function(p) read_snpinfo(panel_dir(p)))

# ---- 1. Structure ---------------------------------------------------------
cat(sprintf("== 1. Structure (chr%d, %s)\n", chr, paste(pops, collapse = ", ")))
ref_info <- infos[[1]]
for (p in pops) {
  info <- infos[[p]]
  h5_snps <- unlist(lapply(panels[[p]], `[[`, "snps"))
  check(identical(h5_snps, info$SNP),
    sprintf("%s: %d blocks list the same %s SNPs, in order, as snpinfo", p, length(panels[[p]]),
      format(nrow(info), big.mark = ",")))
  if (p != pops[1]) {
    check(identical(info[, c("SNP", "BP", "A1", "A2")], ref_info[, c("SNP", "BP", "A1", "A2")]),
      sprintf("%s: same SNPs, positions and alleles as %s", p, pops[1]))
  }
}

# ---- 2. Matrices ----------------------------------------------------------
cat("\n== 2. LD matrices\n")
for (p in pops) {
  stats_by_blk <- t(vapply(panels[[p]], function(b) {
    m <- length(b$snps)
    ld <- b$ld
    if (m == 0) return(c(dim_ok = 1, finite = 1, asym = 0, diag = 0, maxoff = 0, mineig = NA))
    off <- if (m > 1) max(abs(ld[upper.tri(ld)])) else 0
    c(dim_ok = as.numeric(all(dim(ld) == m)),
      finite = as.numeric(all(is.finite(ld))),
      asym = max(abs(ld - t(ld))),
      diag = max(abs(diag(ld) - 1)),
      maxoff = off,
      mineig = min(eigen(ld, symmetric = TRUE, only.values = TRUE)$values))
  }, numeric(6)))
  check(all(stats_by_blk[, "dim_ok"] == 1) && all(stats_by_blk[, "finite"] == 1),
    sprintf("%s: every block is a finite m x m matrix", p))
  check(max(stats_by_blk[, "asym"]) < 1e-8 && max(stats_by_blk[, "diag"]) < 1e-8,
    sprintf("%s: symmetric with unit diagonal (max deviation %.1e)", p,
      max(stats_by_blk[, c("asym", "diag")])))
  check(max(stats_by_blk[, "maxoff"]) <= 1 + 1e-8,
    sprintf("%s: off-diagonal |r| <= 1 (max %.4f)", p, max(stats_by_blk[, "maxoff"])))
  check(min(stats_by_blk[, "mineig"], na.rm = TRUE) > 0,
    sprintf("%s: positive definite (smallest eigenvalue %.2e)", p, min(stats_by_blk[, "mineig"], na.rm = TRUE)))
}

# ---- 3. Agreement with PRS-CS published panels ---------------------------
cat("\n== 3. Agreement with PRS-CS's published 1000 Genomes panels\n")
compare_to_prscs <- function(p, ref_dir) {
  ours <- panels[[p]]
  theirs <- read_panel(ref_dir)
  their_info <- read_snpinfo(ref_dir)
  our_info <- infos[[p]]
  # Allele orientation: +1 same A1/A2, -1 swapped, NA otherwise (dropped).
  m <- match(our_info$SNP, their_info$SNP)
  same <- our_info$A1 == their_info$A1[m] & our_info$A2 == their_info$A2[m]
  swap <- our_info$A1 == their_info$A2[m] & our_info$A2 == their_info$A1[m]
  sgn <- stats::setNames(ifelse(same, 1, ifelse(swap, -1, NA)), our_info$SNP)

  our_blk <- stats::setNames(rep(seq_along(ours), lengths(lapply(ours, `[[`, "snps"))),
    unlist(lapply(ours, `[[`, "snps")))
  our_pos <- stats::setNames(unlist(lapply(ours, function(b) seq_along(b$snps))),
    unlist(lapply(ours, `[[`, "snps")))

  r_ours <- numeric(0)
  r_theirs <- numeric(0)
  for (tb in theirs) {
    keep <- which(tb$snps %in% names(our_blk) & !is.na(sgn[tb$snps]))
    if (length(keep) < 2) next
    s <- tb$snps[keep]
    for (ob in unique(our_blk[s])) {
      k <- keep[our_blk[s] == ob]
      if (length(k) < 2) next
      ss <- tb$snps[k]
      ut <- upper.tri(diag(length(k)))
      a <- ours[[ob]]$ld[our_pos[ss], our_pos[ss], drop = FALSE]
      b <- tb$ld[k, k, drop = FALSE] * outer(sgn[ss], sgn[ss])
      r_ours <- c(r_ours, a[ut])
      r_theirs <- c(r_theirs, b[ut])
    }
  }
  list(n = length(r_ours), cor = stats::cor(r_ours, r_theirs),
    rmse = sqrt(mean((r_ours - r_theirs)^2)), n_snps = sum(!is.na(sgn)))
}
for (p in pops) {
  if (is.na(prscs_root) || !(p %in% names(ref_map))) {
    cat(sprintf("  [SKIP] %s: no published PRS-CS panel to compare with\n", p))
    next
  }
  ref_dir <- file.path(prscs_root, sprintf("ldblk_1kg_%s", ref_map[[p]]))
  res <- compare_to_prscs(p, ref_dir)
  cat(sprintf("  [INFO] %s vs PRS-CS 1KG %s: r correlation %.3f, RMSE %.3f over %s shared within-block pairs (%s SNPs in both)\n",
    p, toupper(ref_map[[p]]), res$cor, res$rmse, format(res$n, big.mark = ","), format(res$n_snps, big.mark = ",")))
}

# ---- 4. Population contrasts ----------------------------------------------
cat("\n== 4. Population contrasts\n")
cat(sprintf("  %-4s %8s %22s %18s %16s\n", "pop", "n", "mean r^2, neighbours", "median MAF", "MAF < 0.05"))
desc <- if (is.na(pop_desc_file)) NULL else
  utils::read.table(pop_desc_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
for (p in pops) {
  r2_nb <- unlist(lapply(panels[[p]], function(b) {
    m <- length(b$snps)
    if (m < 2) return(numeric(0))
    b$ld[cbind(1:(m - 1), 2:m)]^2
  }))
  n <- if (is.null(desc)) NA else sum(desc[[2]][toupper(desc[[3]]) == p])
  maf <- infos[[p]]$MAF
  cat(sprintf("  %-4s %8s %22.3f %18.3f %15.1f%%\n", p, format(n, big.mark = ","), mean(r2_nb),
    stats::median(maf), 100 * mean(maf < 0.05)))
}

cat(sprintf("\nOverall structural and matrix checks: %s\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1)
