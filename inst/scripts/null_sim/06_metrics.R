#!/usr/bin/env Rscript
# Null simulation, step 6: per-replicate, per-arm metrics (plan section 7.1).
#
#   LD fit (no PRS-CS involved), per mixing weight w and scoring
#   regularization delta, over the analysis SNPs of every block b:
#     R_b   = w R_EUR,b + (1 - w) R_ASN,b, then (1 - delta) R_b + delta I
#     L     = sum_b [log det R_b + z_b' R_b^-1 z_b] / m      (primary)
#     C     = sum_b z_b' R_b^-1 z_b / m                     (secondary)
#   For each PRS-CS output (standardized betas, --beta_std=TRUE):
#     score_var_<POP>  variance (denominator n) of Z_T beta over that ancestry's
#                      target individuals, Z_T standardized with discovery constants
#     abs_beta_q50/q90/q99/q999, abs_beta_max, sum_beta2
#     phi              posterior global shrinkage estimate (auto runs only)
#   Identity check: max |beta| difference of mix_w1.0 vs single_EUR and
#   mix_w0.0 vs single_ASN (expected exactly 0).
#
# Usage:
#   Rscript 06_metrics.R <sim_dir> <chr> [pops=EUR,ASN] [deltas=0.001,0]
#
# Reads <sim_dir>/{geno,null,panels/eurcuts,prscs}; writes
# <sim_dir>/metrics/{ld_fit,prscs_metrics,identity}.tsv.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))

args <- ns_args("Rscript 06_metrics.R <sim_dir> <chr> [pops=EUR,ASN] [deltas=0.001,0]", 2)
sim_dir <- args[[1]]
chr <- as.integer(args[[2]])
pops <- toupper(strsplit(ns_arg(args, 3, "EUR,ASN"), ",")[[1]])
deltas <- as.numeric(strsplit(ns_arg(args, 4, "0.001,0"), ",")[[1]])
stopifnot(length(pops) == 2)  # the weight sweep is w on pops[1], 1 - w on pops[2]

geno_dir <- file.path(sim_dir, "geno")
panel_root <- file.path(sim_dir, "panels", "eurcuts")
out_dir <- file.path(sim_dir, "metrics")
dir.create(out_dir, showWarnings = FALSE)

snps <- ns_read_tsv(file.path(geno_dir, "analysis_snps.tsv"))
m <- nrow(snps)

z_files <- sort(list.files(file.path(sim_dir, "null", "z"), "^rep_[0-9]+\\.rds$", full.names = TRUE))
reps <- as.integer(sub("^rep_([0-9]+)\\.rds$", "\\1", basename(z_files)))
zmat <- vapply(z_files, function(f) {
  z <- readRDS(f)
  stopifnot(identical(names(z), snps$SNP))
  z
}, numeric(m))
rownames(zmat) <- snps$SNP
message(sprintf("%d replicates x %d SNPs", ncol(zmat), m))

# ---- LD fit --------------------------------------------------------------
panels <- lapply(stats::setNames(pops, pops), function(p) {
  ns_subset_blocks(ns_read_panel(ns_panel_dir(panel_root, p), chr), snps$SNP)
})
stopifnot(identical(lapply(panels[[1]], `[[`, "snps"), lapply(panels[[2]], `[[`, "snps")))
m_scored <- sum(vapply(panels[[1]], function(b) length(b$snps), integer(1)))
if (m_scored != m) message(sprintf("Note: %d analysis SNPs fall outside every panel block", m - m_scored))

weights <- round(seq(1, 0, by = -0.1), 1)
fit <- list()
for (w in weights) {
  for (delta in deltas) {
    logdet <- 0
    quad <- numeric(ncol(zmat))
    for (b in seq_along(panels[[1]])) {
      s <- panels[[1]][[b]]$snps
      R <- w * panels[[1]][[b]]$ld + (1 - w) * panels[[2]][[b]]$ld
      R <- (1 - delta) * R + delta * diag(length(s))
      U <- chol(R)
      logdet <- logdet + 2 * sum(log(diag(U)))
      v <- backsolve(U, zmat[s, , drop = FALSE], transpose = TRUE)
      quad <- quad + colSums(v^2)
    }
    fit[[length(fit) + 1]] <- data.frame(rep = reps, w = w, delta = delta,
                                         L = (logdet + quad) / m_scored, C = quad / m_scored)
  }
  message(sprintf("LD fit done for w = %.1f", w))
}
fit <- do.call(rbind, fit)
ns_write_tsv(fit, file.path(out_dir, "ld_fit.tsv"))

# ---- PRS-CS outputs --------------------------------------------------------
targets <- lapply(stats::setNames(pops, pops), function(p) ns_load_standardized(geno_dir, p, "target"))

met <- list()
ident <- list()
for (phi_dir in list.dirs(file.path(sim_dir, "prscs"), recursive = FALSE)) {
  phi_mode <- sub("^phi_", "", basename(phi_dir))
  for (rep_dir in list.dirs(phi_dir, recursive = FALSE)) {
    rep <- as.integer(sub("^rep_", "", basename(rep_dir)))
    files <- list.files(rep_dir, sprintf("_pst_eff_.*_chr%d\\.txt$", chr), full.names = TRUE)
    betas <- list()
    for (f in files) {
      arm <- sub("_pst_eff_.*$", "", basename(f))
      beta <- ns_read_beta(f, snps)
      betas[[arm]] <- beta
      ab <- abs(beta)
      row <- data.frame(phi_mode = phi_mode, rep = rep, arm = arm,
                        w = ns_arm_weight(arm, pops),
                        phi = ns_read_phi(file.path(rep_dir, paste0(arm, ".log"))),
                        abs_beta_q50 = stats::quantile(ab, 0.5, names = FALSE),
                        abs_beta_q90 = stats::quantile(ab, 0.9, names = FALSE),
                        abs_beta_q99 = stats::quantile(ab, 0.99, names = FALSE),
                        abs_beta_q999 = stats::quantile(ab, 0.999, names = FALSE),
                        abs_beta_max = max(ab), sum_beta2 = sum(beta^2))
      for (p in pops) {
        sc <- as.vector(targets[[p]] %*% beta)
        row[[paste0("score_var_", p)]] <- mean((sc - mean(sc))^2)
      }
      met[[length(met) + 1]] <- row
    }
    pairs <- list(c("mix_w1.0", paste0("single_", pops[1])), c("mix_w0.0", paste0("single_", pops[2])))
    for (pr in pairs) {
      if (all(pr %in% names(betas))) {
        ident[[length(ident) + 1]] <- data.frame(phi_mode = phi_mode, rep = rep, mix = pr[1], single = pr[2],
                                                 max_abs_diff = max(abs(betas[[pr[1]]] - betas[[pr[2]]])))
      }
    }
  }
}
met <- do.call(rbind, met)
ns_write_tsv(met, file.path(out_dir, "prscs_metrics.tsv"))
ident <- if (length(ident)) do.call(rbind, ident) else data.frame()
ns_write_tsv(ident, file.path(out_dir, "identity.tsv"))
message(sprintf("Wrote %d PRS-CS metric rows; identity check max |diff| = %s",
                nrow(met), if (nrow(ident)) format(max(ident$max_abs_diff)) else "NA"))
