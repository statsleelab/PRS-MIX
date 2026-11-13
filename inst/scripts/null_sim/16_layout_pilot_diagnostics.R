#!/usr/bin/env Rscript
# Layout pilot, sampler-free diagnostics. They bear on two competing
# explanations for the 250 kb accuracy loss and for minority-only LD beating
# the generating mixture:
#
#   Boundary truncation: block-diagonal LD discards correlation between SNPs in
#   different blocks. Measured on the GENERATING covariance of each ancestry,
#   G_k = Z_k'Z_k / n_k over the discovery individuals (exactly the covariance
#   05_simulate_nonnull.R and 04_simulate_null.R draw z from), as the share of
#   squared off-diagonal correlation, over SNP pairs within max_dist, that
#   falls between rather than within blocks. Also given with each pair's r^2
#   reduced by 1/n_k, the expected r^2 between unlinked SNPs in n_k
#   individuals, because at small n_k (e.g. 205 AFR discovery individuals)
#   much of G_k's off-diagonal mass is sampling noise that is nevertheless part
#   of the simulated covariance.
#   Reference mismatch: a panel built from other individuals differs from G_k
#   even within blocks. Measured per block as the squared Frobenius distance
#   between panel and G_k off the diagonal, relative to G_k's off-diagonal mass,
#   for the reference panels and for the discovery-sample panels (which differ
#   from G_k only by the panel estimator), per ancestry and for the mixture at
#   the generating weight.
#
#   LD fit: the Gaussian loss of 06_metrics.R (ridge delta) of each layout x
#   panel source at every weight, on the frozen null z-scores, plus the
#   block-diagonal generating covariance w G_A + (1 - w) G_B ("generating")
#   as the within-block oracle.
#
# Usage:
#   Rscript 16_layout_pilot_diagnostics.R <pilot_dir> <chr> <pops=EUR,AFR> <w0> \
#     <layouts=fixed250kb,fixed1000kb,eurcuts> [max_dist=3000000] [delta=0.001]
#
# Reads <pilot_dir>/{inputs/geno,inputs/null/z,snps/common.tsv,panels/<layout>/<source>};
# writes <pilot_dir>/diagnostics/{blocks,boundary,mismatch,ld_fit}.tsv.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args("Rscript 16_layout_pilot_diagnostics.R <pilot_dir> <chr> <pops> <w0> <layouts> [max_dist] [delta]", 5)
pilot <- args[[1]]
chr <- as.integer(args[[2]])
pops <- toupper(strsplit(args[[3]], ",")[[1]])
w0 <- as.numeric(args[[4]])
layouts <- strsplit(args[[5]], ",")[[1]]
max_dist <- as.numeric(ns_arg(args, 6, "3000000"))
delta <- as.numeric(ns_arg(args, 7, "0.001"))
stopifnot(length(pops) == 2)
sources <- c("reference", "discovery")
out <- file.path(pilot, "diagnostics")
dir.create(out, showWarnings = FALSE)

geno_dir <- file.path(pilot, "inputs", "geno")
all_snps <- ns_read_tsv(file.path(geno_dir, "analysis_snps.tsv"))
snps <- ns_read_tsv(file.path(pilot, "snps", "common.tsv"))
cols <- match(snps$SNP, all_snps$SNP)
stopifnot(!anyNA(cols), !is.unsorted(snps$BP))
m <- nrow(snps)
Z <- lapply(stats::setNames(pops, pops), function(p) ns_load_standardized(geno_dir, p, "discovery")[, cols, drop = FALSE])
n <- vapply(Z, nrow, numeric(1))

# ---- Panels and block membership --------------------------------------------
panels <- list()
block_of <- list()
blocks_tab <- list()
for (layout in layouts) {
  for (src in sources) for (p in pops) {
    blk <- ns_subset_blocks(ns_read_panel(ns_panel_dir(file.path(pilot, "panels", layout, src), p), chr), snps$SNP)
    panels[[layout]][[src]][[p]] <- blk
  }
  ref <- lapply(panels[[layout]][["reference"]][[pops[1]]], `[[`, "snps")
  for (src in sources) for (p in pops) {
    if (!identical(lapply(panels[[layout]][[src]][[p]], `[[`, "snps"), ref)) {
      stop(sprintf("%s panels do not share blocks and SNP order across sources and ancestries", layout), call. = FALSE)
    }
  }
  b <- integer(m)
  for (i in seq_along(ref)) b[match(ref[[i]], snps$SNP)] <- i
  stopifnot(all(b > 0))
  block_of[[layout]] <- b
  sz <- lengths(ref)
  span <- vapply(ref, function(s) diff(range(snps$BP[match(s, snps$SNP)])), numeric(1))
  blocks_tab[[length(blocks_tab) + 1]] <- data.frame(layout = layout, n_blocks = length(ref), m = m,
    snps_median = stats::median(sz), snps_max = max(sz), span_median_bp = stats::median(span), span_max_bp = max(span))
}
ns_write_tsv(do.call(rbind, blocks_tab), file.path(out, "blocks.tsv"))

# ---- Boundary truncation on the generating covariance -----------------------
chunk <- floor(snps$BP / 1e6)
reach <- ceiling(max_dist / 1e6)
bound <- list()
for (p in pops) {
  acc <- matrix(0, length(layouts), 4, dimnames = list(layouts, c("within", "cross", "within_adj", "cross_adj")))
  for (ci in unique(chunk)) {
    I <- which(chunk == ci)
    J <- which(chunk >= ci & chunk <= ci + reach)
    C2 <- (crossprod(Z[[p]][, I, drop = FALSE], Z[[p]][, J, drop = FALSE]) / n[[p]])^2
    ok <- abs(outer(snps$BP[I], snps$BP[J], "-")) <= max_dist & outer(I, J, "<")
    for (layout in layouts) {
      same <- outer(block_of[[layout]][I], block_of[[layout]][J], "==")
      acc[layout, ] <- acc[layout, ] + c(sum(C2[ok & same]), sum(C2[ok & !same]),
                                         sum(C2[ok & same] - 1 / n[[p]]), sum(C2[ok & !same] - 1 / n[[p]]))
    }
  }
  bound[[p]] <- data.frame(pop = p, n_discovery = n[[p]], layout = layouts, max_dist = max_dist, acc,
                           cross_share = acc[, "cross"] / (acc[, "within"] + acc[, "cross"]),
                           cross_share_adj = acc[, "cross_adj"] / (acc[, "within_adj"] + acc[, "cross_adj"]),
                           row.names = NULL)
}
ns_write_tsv(do.call(rbind, bound), file.path(out, "boundary.tsv"))

# ---- Within-block mismatch of panels against the generating covariance ------
gram_block <- function(p, s) {
  k <- match(s, snps$SNP)
  crossprod(Z[[p]][, k, drop = FALSE]) / n[[p]]
}
offdiag2 <- function(x) sum(x^2) - sum(diag(x)^2)
mis <- list()
for (layout in layouts) {
  nb <- length(panels[[layout]][["reference"]][[pops[1]]])
  acc <- list()
  for (b in seq_len(nb)) {
    s <- panels[[layout]][["reference"]][[pops[1]]][[b]]$snps
    G <- lapply(stats::setNames(pops, pops), function(p) gram_block(p, s))
    Gmix <- w0 * G[[1]] + (1 - w0) * G[[2]]
    for (src in sources) {
      R <- lapply(stats::setNames(pops, pops), function(p) panels[[layout]][[src]][[p]][[b]]$ld)
      Rmix <- w0 * R[[1]] + (1 - w0) * R[[2]]
      for (tgt in c(pops, "mix")) {
        key <- paste(src, tgt)
        e <- if (tgt == "mix") offdiag2(Rmix - Gmix) else offdiag2(R[[tgt]] - G[[tgt]])
        d <- if (tgt == "mix") offdiag2(Gmix) else offdiag2(G[[tgt]])
        acc[[key]] <- (if (is.null(acc[[key]])) 0 else acc[[key]]) + c(e, d)
      }
    }
  }
  for (key in names(acc)) {
    kk <- strsplit(key, " ")[[1]]
    mis[[length(mis) + 1]] <- data.frame(layout = layout, source = kk[1], ld = kk[2],
                                         err_offdiag2 = acc[[key]][1], gen_offdiag2 = acc[[key]][2],
                                         rel_err = acc[[key]][1] / acc[[key]][2])
  }
}
ns_write_tsv(do.call(rbind, mis), file.path(out, "mismatch.tsv"))

# ---- LD fit on the frozen null z-scores -------------------------------------
z_files <- sort(list.files(file.path(pilot, "inputs", "null", "z"), "^rep_[0-9]+\\.rds$", full.names = TRUE))
reps <- as.integer(sub("^rep_([0-9]+)\\.rds$", "\\1", basename(z_files)))
zmat <- vapply(z_files, function(f) {
  z <- readRDS(f)
  stopifnot(identical(names(z), all_snps$SNP))
  z[snps$SNP]
}, numeric(m))
rownames(zmat) <- snps$SNP
weights <- round(seq(1, 0, by = -0.1), 1)
fit <- list()
for (layout in layouts) {
  blk <- panels[[layout]][["reference"]][[pops[1]]]
  Gb <- lapply(blk, function(b) lapply(stats::setNames(pops, pops), function(p) gram_block(p, b$snps)))
  for (src in c(sources, "generating")) for (w in weights) {
    logdet <- 0
    quad <- numeric(ncol(zmat))
    for (b in seq_along(blk)) {
      s <- blk[[b]]$snps
      R <- if (src == "generating") w * Gb[[b]][[1]] + (1 - w) * Gb[[b]][[2]]
           else w * panels[[layout]][[src]][[pops[1]]][[b]]$ld + (1 - w) * panels[[layout]][[src]][[pops[2]]][[b]]$ld
      R <- (1 - delta) * R + delta * diag(length(s))
      U <- chol(R)
      logdet <- logdet + 2 * sum(log(diag(U)))
      v <- backsolve(U, zmat[s, , drop = FALSE], transpose = TRUE)
      quad <- quad + colSums(v^2)
    }
    fit[[length(fit) + 1]] <- data.frame(layout = layout, source = src, rep = reps, w = w, delta = delta,
                                         L = (logdet + quad) / m, C = quad / m)
  }
  message(sprintf("LD fit done for %s", layout))
}
ns_write_tsv(do.call(rbind, fit), file.path(out, "ld_fit.tsv"))
message(sprintf("Diagnostics written to %s (%d SNPs, %d null replicates)", out, m, length(reps)))
