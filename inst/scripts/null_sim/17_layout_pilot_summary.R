#!/usr/bin/env Rscript
# Layout pilot summary (see LAYOUT_PILOT.md). Exploratory, not curated paper
# evidence. Reads each pilot setting's metrics, diagnostics, PRS-CS outputs and
# /usr/bin/time records, and writes one page of tables plus every number in TSV.
#
#   1. Inventory: PRS-CS arms found against those expected.
#   2. Reproducibility: the rebuilt 250 kb reference-panel arms against the
#      frozen sweep outputs for the same replicates (posterior effects on the
#      shared SNPs).
#   3. Chain length: long (LONG_ITER) minus default, and xlong minus long, per
#      arm and layout, against the per-replicate SD of the seed rerun; learned phi.
#   4. Contrasts: paired r^2 differences and geometric-mean MSE change by layout,
#      target and panel source (PRS-MIX at w0 against each single panel,
#      discovery against reference panels, the seed rerun).
#   5. Layouts: 1 Mb and EUR-derived cut points minus 250 kb, per arm.
#   6. Weights: mean r^2 at every PRS-MIX weight, per layout.
#   7. Sampler-free diagnostics from 16_layout_pilot_diagnostics.R.
#   8. Run time and peak memory.
# Intervals are pointwise 95% t intervals over replicates, paired within
# replicate, and are not adjusted for multiplicity.
#
# Usage:
#   Rscript 17_layout_pilot_summary.R <pilot_root> <out_dir> [settings=eur_afr_w0.5,eur_asn_w0.5] [chr=22]
#   <pilot_root>/<setting>/ holds runs/<layout>/{metrics/nonnull_metrics.tsv,prscs_nonnull/phi_auto/rep_<r>/},
#   diagnostics/{blocks,boundary,mismatch,ld_fit}.tsv and inputs/existing_prscs_nonnull/rep_<r>/.
# Writes <out_dir>/layout_pilot_summary.{html,tsv}.

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
args <- ns_args("Rscript 17_layout_pilot_summary.R <pilot_root> <out_dir> [settings] [chr=22]", 2)
root <- args[[1]]
out_dir <- args[[2]]
settings <- strsplit(ns_arg(args, 3, "eur_afr_w0.5,eur_asn_w0.5"), ",")[[1]]
chr <- as.integer(ns_arg(args, 4, "22"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
layouts <- c("fixed250kb", "fixed1000kb", "eurcuts")
layout_label <- c(fixed250kb = "250 kb", fixed1000kb = "1 Mb", eurcuts = "EUR cuts")

cells <- list()
put <- function(section, setting, layout, target, item, x, scale = "mean") {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
  q <- if (n > 1) stats::qt(0.975, n - 1) else NA_real_
  lo <- m - q * se
  hi <- m + q * se
  if (scale == "log_pct" && n) {
    lo <- 100 * (exp(lo) - 1)
    hi <- 100 * (exp(hi) - 1)
    se <- 100 * exp(m) * se
    m <- 100 * (exp(m) - 1)
  }
  cells[[length(cells) + 1]] <<- data.frame(section = section, setting = setting, layout = layout, target = target,
                                            item = item, n = n, mean = m, se = se, ci_lo = lo, ci_hi = hi,
                                            sd = if (n > 1 && scale == "mean") stats::sd(x) else NA_real_)
}
per_rep <- function(tab, arm, col) {
  x <- tab[tab$arm == arm, ]
  stats::setNames(x[[col]], x$rep)
}
paired <- function(a, b, f = `-`) {
  k <- intersect(names(a), names(b))
  f(a[k], b[k])
}
read_eff <- function(f) {
  x <- utils::read.table(f, stringsAsFactors = FALSE)
  stats::setNames(x[[6]], x[[2]])
}
parse_time <- function(f) {
  l <- readLines(f)
  el <- sub(".*: ", "", grep("Elapsed \\(wall clock\\)", l, value = TRUE))
  if (!length(el)) return(c(minutes = NA_real_, rss_mb = NA_real_))  # arm still running
  p <- as.numeric(strsplit(el, ":")[[1]])
  sec <- sum(p * 60^(rev(seq_along(p)) - 1))
  rss <- as.numeric(sub(".*: ", "", grep("Maximum resident set size", l, value = TRUE)))
  c(minutes = sec / 60, rss_mb = rss / 1024)
}
chain_of <- function(arm) ifelse(grepl("^xlong_", arm), "8000", ifelse(grepl("^long_", arm), "4000", "1000"))

info <- list()
for (s in settings) {
  pair <- sub("_w[0-9.]+$", "", s)
  w0 <- as.numeric(sub("^.*_w", "", s))
  pops <- toupper(strsplit(pair, "_")[[1]])
  A <- pops[1]
  B <- pops[2]
  P <- file.path(root, s)
  mix0 <- sprintf("mix_w%.1f", w0)
  main_arms <- c(paste0("single_", pops), sprintf("mix_w%.1f", seq(0.9, 0.1, by = -0.1)), paste0(mix0, "_seed2"),
                 paste0("single_", A, "_seed2"), paste0("disc_single_", pops), paste0("disc_", mix0))
  long_arms <- c(paste0("single_", pops), mix0)
  met <- do.call(rbind, lapply(layouts, function(l) {
    x <- ns_read_tsv(file.path(P, "runs", l, "metrics", "nonnull_metrics.tsv"))
    x$layout <- l
    x[x$phi_mode == "auto", ]
  }))
  reps <- sort(unique(met$rep[met$arm == paste0("single_", A)]))
  long_reps <- sort(unique(met$rep[met$arm == paste0("long_", mix0)]))
  has_xlong <- any(grepl("^xlong_", met$arm))
  info[[s]] <- list(pair = pair, w0 = w0, A = A, B = B, reps = reps, long_reps = long_reps, has_xlong = has_xlong,
                    xlong_reps = sort(unique(met$rep[met$arm == paste0("xlong_", mix0)])))

  # 1. Inventory
  for (l in layouts) {
    m <- met[met$layout == l, ]
    expect <- c(length(main_arms) * length(reps), 3 * length(long_reps), if (has_xlong) 3 * length(info[[s]]$xlong_reps) else 0)
    found <- c(sum(m$arm %in% main_arms & m$rep %in% reps), sum(m$arm %in% paste0("long_", long_arms)),
               sum(m$arm %in% paste0("xlong_", long_arms)))
    put("inventory", s, l, "", "main_found", found[1]); put("inventory", s, l, "", "main_expected", expect[1])
    put("inventory", s, l, "", "long_found", found[2]); put("inventory", s, l, "", "long_expected", expect[2])
    put("inventory", s, l, "", "xlong_found", found[3]); put("inventory", s, l, "", "xlong_expected", expect[3])
  }

  # 2. Reproducibility of the 250 kb reference arms against the frozen sweep outputs
  run_dir <- file.path(P, "runs", "fixed250kb", "prscs_nonnull", "phi_auto")
  old_dir <- file.path(P, "inputs", "existing_prscs_nonnull")
  for (arm in c(paste0("single_", pops), mix0, paste0(mix0, "_seed2"), paste0("single_", A, "_seed2"))) {
    stat <- t(vapply(reps, function(r) {
      tail <- sprintf("%s_pst_eff_a1_b0.5_phiauto_chr%d.txt", arm, chr)
      fn <- file.path(run_dir, sprintf("rep_%03d", r), tail)
      fo <- file.path(old_dir, sprintf("rep_%03d", r), tail)
      if (!file.exists(fn) || !file.exists(fo)) return(c(NA, NA, NA))
      a <- read_eff(fn)
      b <- read_eff(fo)
      k <- intersect(names(a), names(b))
      c(max(abs(a[k] - b[k])) / max(abs(b[k])), stats::cor(a[k], b[k]), length(union(names(a), names(b))) - length(k))
    }, numeric(3)))
    put("repro", s, "fixed250kb", "", paste0(arm, ":max_rel_diff"), stat[, 1])
    put("repro", s, "fixed250kb", "", paste0(arm, ":cor"), stat[, 2])
    put("repro", s, "fixed250kb", "", paste0(arm, ":snps_not_shared"), stat[, 3])
  }

  for (l in layouts) {
    m <- met[met$layout == l, ]
    # 3. Chain length
    for (tg in pops) {
      col <- paste0("r2_", tg)
      seed_sd <- stats::sd(paired(per_rep(m, paste0(mix0, "_seed2"), col), per_rep(m, mix0, col)))
      put("chain", s, l, tg, "seed_sd_per_rep", seed_sd)
      for (arm in long_arms) {
        def <- per_rep(m, arm, col)
        lng <- per_rep(m, paste0("long_", arm), col)
        put("chain", s, l, tg, paste0(arm, ":long_minus_default"), paired(lng, def))
        if (has_xlong) put("chain", s, l, tg, paste0(arm, ":xlong_minus_long"), paired(per_rep(m, paste0("xlong_", arm), col), lng))
      }
      # Does chain length move the PRS-MIX contrasts?
      for (ch in c("", "long_", if (has_xlong) "xlong_")) {
        r <- if (ch == "") long_reps else info[[s]][[if (ch == "long_") "long_reps" else "xlong_reps"]]
        for (sg in pops) {
          d <- paired(per_rep(m, paste0(ch, mix0), col), per_rep(m, paste0(ch, "single_", sg), col))
          put("chain_contrast", s, l, tg, sprintf("%s:mix_minus_single_%s", chain_of(paste0(ch, "x")), sg), d[names(d) %in% r])
        }
      }
    }
    for (arm in long_arms) {
      put("phi", s, l, "", paste0(arm, ":1000"), per_rep(m, arm, "phi")[as.character(long_reps)])
      put("phi", s, l, "", paste0(arm, ":4000"), per_rep(m, paste0("long_", arm), "phi"))
      if (has_xlong) put("phi", s, l, "", paste0(arm, ":8000"), per_rep(m, paste0("xlong_", arm), "phi"))
    }

    # 4. Contrasts at w0
    for (tg in pops) {
      r2 <- paste0("r2_", tg)
      ms <- paste0("mse_gv_", tg)
      cmp <- list(c(mix0, paste0("single_", A)), c(mix0, paste0("single_", B)),
                  c(paste0("disc_", mix0), paste0("disc_single_", A)), c(paste0("disc_", mix0), paste0("disc_single_", B)),
                  c(paste0("disc_", mix0), mix0), c(paste0("disc_single_", A), paste0("single_", A)),
                  c(paste0("disc_single_", B), paste0("single_", B)),
                  c(paste0(mix0, "_seed2"), mix0), c(paste0("single_", A, "_seed2"), paste0("single_", A)))
      for (cc in cmp) {
        nm <- paste(cc, collapse = " - ")
        put("contrast_r2", s, l, tg, nm, paired(per_rep(m, cc[1], r2), per_rep(m, cc[2], r2)))
        put("contrast_mse", s, l, tg, nm, paired(per_rep(m, cc[1], ms), per_rep(m, cc[2], ms), function(a, b) log(a / b)), "log_pct")
      }
      for (arm in c(paste0("single_", pops), mix0, paste0("disc_single_", pops), paste0("disc_", mix0))) {
        put("r2", s, l, tg, arm, per_rep(m, arm, r2))
      }
      # 6. Weights
      for (w in c(1, seq(0.9, 0.1, by = -0.1), 0)) {
        arm <- if (w == 1) paste0("single_", A) else if (w == 0) paste0("single_", B) else sprintf("mix_w%.1f", w)
        put("weight_r2", s, l, tg, sprintf("%.1f", w), per_rep(m, arm, r2))
      }
    }
  }

  # 5. Layouts against 250 kb
  base <- met[met$layout == "fixed250kb", ]
  for (l in setdiff(layouts, "fixed250kb")) {
    m <- met[met$layout == l, ]
    for (tg in pops) for (arm in c(paste0("single_", pops), mix0, paste0("disc_single_", pops), paste0("disc_", mix0))) {
      col <- paste0("r2_", tg)
      put("layout_r2", s, l, tg, arm, paired(per_rep(m, arm, col), per_rep(base, arm, col)))
    }
  }

  # 7. Diagnostics
  dg <- file.path(P, "diagnostics")
  bd <- ns_read_tsv(file.path(dg, "boundary.tsv"))
  for (i in seq_len(nrow(bd))) {
    put("boundary", s, bd$layout[i], bd$pop[i], "cross_share", bd$cross_share[i])
    put("boundary", s, bd$layout[i], bd$pop[i], "cross_share_adj", bd$cross_share_adj[i])
  }
  mi <- ns_read_tsv(file.path(dg, "mismatch.tsv"))
  for (i in seq_len(nrow(mi))) put("mismatch", s, mi$layout[i], mi$ld[i], mi$source[i], mi$rel_err[i])
  fit <- ns_read_tsv(file.path(dg, "ld_fit.tsv"))
  for (l in layouts) for (src in unique(fit$source)) {
    f <- fit[fit$layout == l & fit$source == src, ]
    best <- vapply(split(f, f$rep), function(x) x$w[which.min(x$L)], numeric(1))
    put("ld_fit_argmin", s, l, "", src, best)
  }

  # 8. Run time and memory
  for (l in layouts) {
    tf <- Sys.glob(file.path(P, "runs", l, "prscs_nonnull", "phi_auto", "rep_*", "*.time"))
    if (!length(tf)) next
    tm <- t(vapply(tf, parse_time, numeric(2)))
    ok <- is.finite(tm[, "minutes"])
    tm <- tm[ok, , drop = FALSE]
    ch <- chain_of(sub("\\.time$", "", basename(tf[ok])))
    for (k in unique(ch)) {
      put("runtime", s, l, k, "minutes", tm[ch == k, "minutes"])
      put("runtime", s, l, k, "peak_rss_mb", tm[ch == k, "rss_mb"])
    }
  }
}
cells <- do.call(rbind, cells)
ns_write_tsv(cells, file.path(out_dir, "layout_pilot_summary.tsv"))

# ---- Page ------------------------------------------------------------------
get1 <- function(section, s, l, tg, item) {
  r <- cells[cells$section == section & cells$setting == s & cells$layout == l & cells$target == tg & cells$item == item, ]
  if (nrow(r)) r[1, ] else NULL
}
num <- function(x, d, signed = FALSE) {
  out <- sub("^-", "&minus;", formatC(x, format = "f", digits = d, flag = if (signed) "+" else ""))
  ifelse(is.na(x), "&ndash;", out)
}
ci <- function(r, d = 4, signed = TRUE) {
  if (is.null(r) || r$n == 0) return("&ndash;")
  if (r$n == 1) return(num(r$mean, d, signed))
  sprintf("%s <span class=\"ci\">[%s, %s]</span>", num(r$mean, d, signed), num(r$ci_lo, d, signed), num(r$ci_hi, d, signed))
}
tab <- function(head, rows) {
  c("<div class=\"wrap\"><table>", "<thead><tr>", sprintf("<th>%s</th>", head), "</tr></thead>", "<tbody>",
    vapply(rows, function(r) paste0("<tr>", paste(sprintf("<td>%s</td>", r), collapse = ""), "</tr>"), ""),
    "</tbody></table></div>")
}
sec <- function(title, caption, body) c("<section>", sprintf("<h2>%s</h2>", title), sprintf("<p class=\"cap\">%s</p>", caption), body, "</section>")
by_layout <- function(f) unname(vapply(layouts, f, ""))
lh <- unname(layout_label[layouts])
pretty_arm <- function(arm, I) {
  arm <- sub("_seed2$", " (seed 2)", arm)
  arm <- sub(sprintf("^disc_mix_w%.1f", I$w0), "PRS-MIX w<sub>0</sub>, discovery", arm)
  arm <- sub(sprintf("^mix_w%.1f", I$w0), "PRS-MIX w<sub>0</sub>", arm)
  arm <- sub("^disc_single_(.*)$", "\\1-only, discovery", arm)
  sub("^single_(.*)$", "\\1-only", arm)
}

body <- character(0)
for (s in settings) {
  I <- info[[s]]
  pops <- c(I$A, I$B)
  mix0 <- sprintf("mix_w%.1f", I$w0)
  title <- sprintf("%s+%s %g:%g", I$A, I$B, 100 * I$w0, 100 * (1 - I$w0))
  body <- c(body, sprintf("<h1>%s</h1>", title))

  inv <- by_layout(function(l) {
    g <- function(k) get1("inventory", s, l, "", k)$mean
    sprintf("%d/%d main, %d/%d long, %d/%d 8,000-iteration", g("main_found"), g("main_expected"), g("long_found"),
            g("long_expected"), g("xlong_found"), g("xlong_expected"))
  })
  rep_rows <- lapply(c(paste0("single_", pops), mix0, paste0(mix0, "_seed2"), paste0("single_", I$A, "_seed2")), function(arm) {
    a <- get1("repro", s, "fixed250kb", "", paste0(arm, ":max_rel_diff"))
    b <- get1("repro", s, "fixed250kb", "", paste0(arm, ":cor"))
    k <- get1("repro", s, "fixed250kb", "", paste0(arm, ":snps_not_shared"))
    c(pretty_arm(arm, I), a$n, formatC(a$mean, format = "g", digits = 2),
      formatC(b$mean, format = "f", digits = 6), num(k$mean, 1))
  })
  body <- c(body, sec(sprintf("%s. Inventory and reproducibility", title),
    paste("Arms found by layout:", paste(sprintf("%s: %s", lh, inv), collapse = "; "),
          ". Reproducibility: the pilot's rebuilt 250 kb reference-panel arms against the frozen sweep outputs for the same replicates, on the SNPs both share; largest absolute effect difference relative to the largest sweep effect, and correlation, averaged over replicates."),
    tab(c("Arm", "Replicates", "Max relative difference", "Correlation", "SNPs not shared"), rep_rows)))

  # Chain length
  ch_rows <- list()
  for (tg in pops) for (arm in c(paste0("single_", pops), mix0)) {
    ch_rows[[length(ch_rows) + 1]] <- c(tg, pretty_arm(arm, I), by_layout(function(l) {
      a <- get1("chain", s, l, tg, paste0(arm, ":long_minus_default"))
      b <- get1("chain", s, l, tg, paste0(arm, ":xlong_minus_long"))
      paste0(ci(a), if (!is.null(b) && b$n) paste0("<br>", ci(b)) else "")
    }))
  }
  sd_row <- by_layout(function(l) paste(vapply(pops, function(tg) num(get1("chain", s, l, tg, "seed_sd_per_rep")$mean, 4), ""), collapse = " / "))
  phi_rows <- lapply(c(paste0("single_", pops), mix0), function(arm) c(pretty_arm(arm, I), by_layout(function(l) {
    v <- vapply(c("1000", "4000", "8000"), function(k) {
      r <- get1("phi", s, l, "", paste0(arm, ":", k))
      if (is.null(r) || r$n == 0) "&ndash;" else formatC(r$mean, format = "g", digits = 3)
    }, "")
    paste(v, collapse = " / ")
  })))
  cc_rows <- list()
  for (tg in pops) for (sg in pops) {
    cc_rows[[length(cc_rows) + 1]] <- c(tg, sprintf("PRS-MIX w<sub>0</sub> &minus; %s-only", sg), by_layout(function(l) {
      paste(vapply(c("1000", "4000", "8000"), function(k) {
        r <- get1("chain_contrast", s, l, tg, sprintf("%s:mix_minus_single_%s", k, sg))
        if (is.null(r) || r$n == 0) "&ndash;" else num(r$mean, 4, TRUE)
      }, ""), collapse = " / ")
    }))
  }
  body <- c(body, sec(sprintf("%s. Chain length", title),
    sprintf("Paired r<sup>2</sup> change, first line 4,000 minus 1,000 iterations (replicates %s), second line 8,000 minus 4,000 (replicates %s), seed 1000 + rep throughout. Prespecified pass: within the per-replicate SD of the seed rerun (PRS-MIX w<sub>0</sub>, seed 2000 minus 1000; first / second-ancestry targets): %s.",
            paste(range(I$long_reps), collapse = "&ndash;"), if (I$has_xlong) paste(range(I$xlong_reps), collapse = "&ndash;") else "none yet",
            paste(sprintf("%s %s", lh, sd_row), collapse = "; ")),
    c(tab(c("Targets", "Arm", lh), ch_rows),
      "<p class=\"cap\">Learned &phi;, mean over the same replicates, at 1,000 / 4,000 / 8,000 iterations.</p>",
      tab(c("Arm", lh), phi_rows),
      "<p class=\"cap\">Do longer chains move the comparisons? Mean paired r<sup>2</sup> difference on the chain-check replicates at 1,000 / 4,000 / 8,000 iterations.</p>",
      tab(c("Targets", "Contrast", lh), cc_rows))))

  # Contrasts
  cmp <- c(paste(mix0, "-", paste0("single_", pops)), paste(paste0("disc_", mix0), "-", paste0("disc_single_", pops)),
           paste(paste0("disc_", mix0), "-", mix0), paste(paste0("disc_single_", pops), "-", paste0("single_", pops)),
           paste(paste0(mix0, "_seed2"), "-", mix0), paste(paste0("single_", I$A, "_seed2"), "-", paste0("single_", I$A)))
  ct_rows <- list()
  for (tg in pops) for (nm in cmp) {
    parts <- strsplit(nm, " - ")[[1]]
    ct_rows[[length(ct_rows) + 1]] <- c(tg, paste(pretty_arm(parts[1], I), "&minus;", pretty_arm(parts[2], I)),
      by_layout(function(l) paste0(ci(get1("contrast_r2", s, l, tg, nm)), "<br>", ci(get1("contrast_mse", s, l, tg, nm), 1), "%")))
  }
  body <- c(body, sec(sprintf("%s. Paired contrasts", title),
    sprintf("Replicates %d&ndash;%d, default sampler. First line: r<sup>2</sup> difference with the true genetic value (positive favours the first arm). Second line: geometric-mean MSE change (negative favours the first arm). Discovery panels are built from the GWAS's own discovery individuals; they are a diagnostic, not a usable method.",
            min(I$reps), max(I$reps)),
    tab(c("Targets", "Contrast", lh), ct_rows)))

  # Layouts, absolute r^2
  arms <- c(paste0("single_", pops), mix0, paste0("disc_single_", pops), paste0("disc_", mix0))
  ly_rows <- list()
  for (tg in pops) for (arm in arms) {
    ly_rows[[length(ly_rows) + 1]] <- c(tg, pretty_arm(arm, I),
      num(get1("r2", s, "fixed250kb", tg, arm)$mean, 3),
      ci(get1("layout_r2", s, "fixed1000kb", tg, arm)), ci(get1("layout_r2", s, "eurcuts", tg, arm)))
  }
  body <- c(body, sec(sprintf("%s. Layouts", title),
    "Mean r<sup>2</sup> at 250 kb, and the paired change at 1 Mb and at EUR-derived cut points relative to 250 kb, default sampler. 1 Mb tiles and EUR cuts differ in both block size and boundary placement.",
    tab(c("Targets", "Arm", "250 kb r<sup>2</sup>", "1 Mb &minus; 250 kb", "EUR cuts &minus; 250 kb"), ly_rows)))

  # Weights
  ws <- sprintf("%.1f", c(1, seq(0.9, 0.1, by = -0.1), 0))
  w_rows <- list()
  for (tg in pops) for (l in layouts) {
    v <- vapply(ws, function(w) get1("weight_r2", s, l, tg, w)$mean, numeric(1))
    w_rows[[length(w_rows) + 1]] <- c(tg, layout_label[[l]], ifelse(v == max(v), sprintf("<b>%s</b>", num(v, 3)), num(v, 3)))
  }
  body <- c(body, sec(sprintf("%s. Weights", title),
    sprintf("Mean r<sup>2</sup> of reference-panel PRS-MIX at each weight on %s (w = 1 is %s-only, w = 0 is %s-only), default sampler; bold: highest in the row. The best weight is chosen with the true genetic values, so it is an oracle.", I$A, I$A, I$B),
    tab(c("Targets", "Layout", ws), w_rows)))

  # Diagnostics
  bd_rows <- lapply(pops, function(p) c(p, by_layout(function(l) {
    a <- get1("boundary", s, l, p, "cross_share")$mean
    b <- get1("boundary", s, l, p, "cross_share_adj")$mean
    sprintf("%s%% (%s%%)", num(100 * a, 0), num(100 * b, 0))
  })))
  mi_rows <- list()
  for (tg in c(pops, "mix")) for (src in c("reference", "discovery")) {
    mi_rows[[length(mi_rows) + 1]] <- c(if (tg == "mix") "Mixture at w<sub>0</sub>" else tg, src,
      by_layout(function(l) sprintf("%s%%", num(100 * get1("mismatch", s, l, tg, src)$mean, 1))))
  }
  fit_rows <- lapply(c("reference", "discovery", "generating"), function(src) c(src, by_layout(function(l) {
    r <- get1("ld_fit_argmin", s, l, "", src)
    sprintf("%s (%s)", num(r$mean, 2), num(r$se, 2))
  })))
  body <- c(body, sec(sprintf("%s. Sampler-free diagnostics", title),
    "Boundary: share of squared correlation in the generating covariance, over SNP pairs within 3 Mb, that falls between blocks (in brackets, with each pair's r<sup>2</sup> reduced by 1/n of the discovery sample). Mismatch: within-block off-diagonal squared error of each panel against the generating covariance, relative to the generating covariance's off-diagonal mass. LD fit: weight minimising the null-z loss per replicate, mean (SE); generating is the block-diagonal generating covariance.",
    c(tab(c("Ancestry", lh), bd_rows), tab(c("LD", "Panel source", lh), mi_rows), tab(c("Panel source", lh), fit_rows))))

  rt_rows <- list()
  for (k in c("1000", "4000", "8000")) {
    r <- lapply(layouts, function(l) list(t = get1("runtime", s, l, k, "minutes"), m = get1("runtime", s, l, k, "peak_rss_mb")))
    if (all(vapply(r, function(x) is.null(x$t), TRUE))) next
    rt_rows[[length(rt_rows) + 1]] <- c(format(as.integer(k), big.mark = ","), vapply(r, function(x)
      if (is.null(x$t)) "&ndash;" else sprintf("%s min, %s MB (%d arms)", num(x$t$mean, 1), num(x$m$mean, 0), x$t$n), ""))
  }
  body <- c(body, sec(sprintf("%s. Run time and memory", title),
    "Mean wall time and peak resident memory per PRS-CS arm (one CPU), from /usr/bin/time.",
    tab(c("Iterations", lh), rt_rows)))
}

git_head <- suppressWarnings(tryCatch(system2("git", c(paste0("--git-dir=", file.path(script_dir, "..", "..", "..", ".git")), "rev-parse", "HEAD"),
                                             stdout = TRUE, stderr = FALSE), error = function(e) character(0)))
git_head <- if (length(git_head)) substr(git_head[1], 1, 7) else "unknown"
html <- c(
  "<!doctype html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">", "<title>Layout pilot summary</title>",
  "<!-- Generated by inst/scripts/null_sim/17_layout_pilot_summary.R; do not edit. -->", "<style>",
  "body { margin: 0; padding: 40px 24px 56px; background: #fff; color: #111; font: 15px/1.45 Georgia, \"Times New Roman\", Times, serif; }",
  "main { max-width: 1180px; margin: 0 auto; display: flex; flex-direction: column; gap: 36px; }",
  "h1 { font: 700 21px/1.3 -apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif; margin: 16px 0 0; border-top: 2px solid #111; padding-top: 14px; }",
  "h2 { font: 600 17px/1.3 -apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif; margin: 0 0 8px; }",
  ".cap { margin: 0 0 12px; font-size: 14px; }", ".wrap { overflow-x: auto; margin-bottom: 14px; }", ".ci { color: #555; font-size: 12.5px; }",
  "table { border-collapse: collapse; width: 100%; font-variant-numeric: lining-nums tabular-nums; font-size: 14px; border-top: 2px solid #111; border-bottom: 2px solid #111; }",
  "th, td { padding: 4px 10px; text-align: center; white-space: nowrap; vertical-align: top; }",
  "th:nth-child(-n+2), td:nth-child(-n+2) { text-align: left; }",
  "thead th { font-weight: 700; border-bottom: 1px solid #111; }",
  "footer { font: 12px/1.5 -apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif; color: #666; }",
  "</style>", "</head>", "<body>", "<main>",
  "<p class=\"cap\">Exploratory, not curated paper evidence (see <code>inst/scripts/null_sim/LAYOUT_PILOT.md</code>). Chromosome 22, PRS-CS with learned &phi;, a = 1, b = 0.5; default sampler 1,000 iterations, 500 burn-in, thinning 5, seed 1000 + rep unless stated. Square brackets are pointwise 95% t intervals over replicates, paired within replicate, not adjusted for multiplicity.</p>",
  body,
  sprintf("<footer>Generated by <code>inst/scripts/null_sim/17_layout_pilot_summary.R</code> at PRS-MIX %s. Every number is in <code>layout_pilot_summary.tsv</code>.</footer>",
          git_head),
  "</main>", "</body>", "</html>")
writeLines(html, file.path(out_dir, "layout_pilot_summary.html"))
message(sprintf("Wrote %s: %d settings, %d cells", out_dir, length(settings), nrow(cells)))
