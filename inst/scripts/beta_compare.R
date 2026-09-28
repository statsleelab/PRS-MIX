args <- commandArgs(trailingOnly = TRUE)
eur_path <- args[[1]]
blend_paths <- strsplit(args[[2]], ",")[[1]]
blend_labels <- strsplit(args[[3]], ",")[[1]]
out_dir <- args[[4]]
n_sample <- if (length(args) >= 5) as.integer(args[[5]]) else 15000

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_pst <- function(path) {
  d <- utils::read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                         col.names = c("chr", "snp", "bp", "a1", "a2", "beta"))
  d[!duplicated(d$snp), ]
}

eur <- read_pst(eur_path)

summary_rows <- list()

for (i in seq_along(blend_paths)) {
  label <- blend_labels[i]
  blend <- read_pst(blend_paths[i])

  m <- merge(eur, blend, by = "snp", suffixes = c("_eur", "_blend"))
  same_orient <- m$a1_eur == m$a1_blend & m$a2_eur == m$a2_blend
  swapped     <- m$a1_eur == m$a2_blend & m$a2_eur == m$a1_blend
  keep <- same_orient | swapped
  m <- m[keep, ]
  same_orient <- same_orient[keep]
  m$beta_blend_aligned <- ifelse(same_orient, m$beta_blend, -m$beta_blend)

  n_total <- nrow(m)
  fit <- stats::lm(beta_blend_aligned ~ beta_eur, data = m)
  r <- stats::cor(m$beta_eur, m$beta_blend_aligned)

  set.seed(42)
  samp <- m[sample.int(n_total, min(n_sample, n_total)), c("beta_eur", "beta_blend_aligned")]
  names(samp) <- c("beta_eur", "beta_blend")
  utils::write.csv(samp, file.path(out_dir, sprintf("sample_%s.csv", label)), row.names = FALSE)

  summary_rows[[label]] <- data.frame(
    label = label,
    n_snps_matched = n_total,
    r = round(r, 4),
    slope = round(unname(stats::coef(fit)[2]), 4),
    intercept = round(unname(stats::coef(fit)[1]), 6),
    stringsAsFactors = FALSE
  )
  message(sprintf("%s: n=%d matched, r=%.4f, slope=%.4f", label, n_total, r, unname(stats::coef(fit)[2])))
}

summary_df <- do.call(rbind, summary_rows)
utils::write.csv(summary_df, file.path(out_dir, "beta_compare_summary.csv"), row.names = FALSE)
print(summary_df)
