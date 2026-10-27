test_that("computeSuperPopLDWindow validates superpopulation label", {
  td <- tempdir()
  pop_desc <- file.path(td, "pop_desc_sup.txt")

  writeLines(c(
    "Population_Abbreviation\tNumber_of_Subjects\tSuper_Population",
    "P1\t10\tEUR",
    "P2\t10\tAFR"
  ), pop_desc)

  expect_error(
    computeSuperPopLDWindow(
      chr = 1,
      start_bp = 1,
      end_bp = 100,
      superpopulation_label = "XXX",
      reference_index_file = "missing_index.gz",
      reference_data_file = "missing_geno.gz",
      reference_pop_desc_file = pop_desc
    ),
    "invalid population name|Invalid superpopulation"
  )
})

.covariance_fixture <- function() {
  root <- testthat::test_path("fixtures", "ld-covariance")
  list(
    reference_index_file = file.path(root, "index.bgz"),
    reference_data_file = file.path(root, "genotypes.bgz"),
    reference_pop_desc_file = file.path(root, "populations.tsv")
  )
}

.covariance_genotypes <- function() {
  # Rows are individuals, columns are SNPs. These exactly match the BGZF fixture.
  list(
    cbind(c(0, 0, 2, 2), c(0, 0, 2, 2), c(0, 0, 0, 2)),
    cbind(rep(c(0, 0, 2, 2), 3), rep(c(2, 2, 0, 0), 3),
          rep(c(1, 1, 2, 2), 3))
  )
}

.expected_superpop_covariance <- function(genotypes) {
  weights <- vapply(genotypes, nrow, integer(1))
  weights <- weights / sum(weights)
  means <- lapply(genotypes, colMeans)
  overall_mean <- Reduce(`+`, Map(function(w, mu) w * mu, weights, means))
  # Law of total covariance, with sample covariance within each population.
  Reduce(`+`, Map(function(w, g, mu) {
    w * (stats::cov(g) + tcrossprod(mu - overall_mean))
  }, weights, genotypes, means))
}

test_that("superpopulation LD uses normalized sample covariances for unequal populations", {
  genotypes <- .covariance_genotypes()
  expected <- stats::cov2cor(.expected_superpop_covariance(genotypes))
  result <- do.call(computeSuperPopLDWindow, c(
    list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR"),
    .covariance_fixture()
  ))

  expect_identical(result$snplist$rsid, c("rs1", "rs2", "rs3"))
  expect_equal(unname(result$cormat), unname(expected), tolerance = 1e-12)
  # Opposing LD in populations of sizes 4 and 12. The former bug gave -0.913386.
  expect_equal(result$cormat[1, 2], -8 / 19, tolerance = 1e-12)
  # SNP3 has different population means, exercising the between-population term.
  expect_equal(result$snplist$af1pop,
               colMeans(do.call(rbind, genotypes)) / 2, ignore_attr = TRUE)
  expect_gte(min(eigen(result$cormat, symmetric = TRUE)$values), -1e-12)
})

test_that("single-population LD still equals genotype Pearson correlation", {
  genotypes <- .covariance_genotypes()
  for (i in seq_along(genotypes)) {
    result <- do.call(computeSuperPopLDWindow, c(
      list(chr = 1, start_bp = 1, end_bp = 100,
           superpopulation_label = paste0("P", i)),
      .covariance_fixture()
    ))
    expect_equal(unname(result$cormat), unname(stats::cor(genotypes[[i]])),
                 tolerance = 1e-12)
  }
})

test_that("panel writer stores the corrected LD from the compiled reader", {
  skip_if_not_installed("hdf5r")
  panel_dir <- tempfile("covariance-panel-")
  dir.create(panel_dir)
  do.call(write_prscs_ld_chr, c(
    list(chr = 1, superpopulation_label = "EUR", panel_dir = panel_dir,
         block_size = 100, shrinkage = 0, verbose = FALSE),
    .covariance_fixture()
  ))
  h5 <- hdf5r::H5File$new(file.path(panel_dir, "ldblk_1kg_chr1.hdf5"), "r")
  on.exit(h5$close_all(), add = TRUE)
  expected <- stats::cov2cor(.expected_superpop_covariance(.covariance_genotypes()))
  expect_identical(h5[["blk_1/snplist"]]$read(), c("rs1", "rs2", "rs3"))
  expect_equal(unname(h5[["blk_1/ldblk"]]$read()), unname(expected),
               tolerance = 1e-10)
})

test_that("keep_bp restricts the compiled reader to the requested positions", {
  fx <- .covariance_fixture()
  full <- do.call(computeSuperPopLDWindow, c(
    list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR"), fx))
  kept <- do.call(computeSuperPopLDWindow, c(
    list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR",
         keep_bp = c(10L, 30L)), fx))

  expect_identical(kept$snplist$rsid, c("rs1", "rs3"))
  expect_equal(unname(kept$cormat), unname(full$cormat[c(1, 3), c(1, 3)]), tolerance = 1e-12)
  expect_error(
    do.call(computeSuperPopLDWindow, c(
      list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR",
           keep_bp = 999L), fx)),
    "No SNPs found in requested region"
  )
})

test_that("min_snps = 1 returns a single SNP as a 1x1 block", {
  fx <- .covariance_fixture()
  one <- do.call(computeSuperPopLDWindow, c(
    list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR",
         keep_bp = 20L, min_snps = 1L), fx))
  expect_identical(one$snplist$rsid, "rs2")
  expect_equal(unname(one$cormat), matrix(1, 1, 1))
  expect_error(
    do.call(computeSuperPopLDWindow, c(
      list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR",
           keep_bp = 20L), fx)),
    "Fewer than 2 SNPs remain"
  )
})

# sample_keep / readSuperPopGenotypes -------------------------------------

test_that("readSuperPopGenotypes returns the fixture's dosages in population order", {
  genotypes <- .covariance_genotypes()
  result <- do.call(readSuperPopGenotypes, c(
    list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR"),
    .covariance_fixture()
  ))

  expect_identical(result$snplist$rsid, c("rs1", "rs2", "rs3"))
  expect_equal(unname(result$genotypes), unname(do.call(rbind, genotypes)), ignore_attr = TRUE)
  expect_identical(result$samples$population, rep(c("P1", "P2"), c(4, 12)))
  expect_identical(result$samples$position, c(1:4, 1:12))
})

test_that("sample_keep restricts genotypes and LD to the kept individuals", {
  genotypes <- .covariance_genotypes()
  keep <- list(P1 = c(1L, 3L, 4L), P2 = c(2L, 5L, 6L, 9L, 12L))
  kept <- list(genotypes[[1]][keep$P1, , drop = FALSE], genotypes[[2]][keep$P2, , drop = FALSE])

  geno <- do.call(readSuperPopGenotypes, c(
    list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR", sample_keep = keep),
    .covariance_fixture()
  ))
  expect_equal(unname(geno$genotypes), unname(do.call(rbind, kept)), ignore_attr = TRUE)
  expect_identical(geno$samples$position, c(keep$P1, keep$P2))

  ld <- do.call(computeSuperPopLDWindow, c(
    list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR", sample_keep = keep),
    .covariance_fixture()
  ))
  # Population weights follow the kept counts (3 and 5), not the full sizes.
  expect_equal(unname(ld$cormat), unname(stats::cov2cor(.expected_superpop_covariance(kept))),
               tolerance = 1e-12)
  expect_equal(ld$snplist$af1pop, colMeans(do.call(rbind, kept)) / 2, ignore_attr = TRUE)
})

test_that("sample_keep with every individual matches no sample_keep", {
  full <- do.call(computeSuperPopLDWindow, c(
    list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR"),
    .covariance_fixture()
  ))
  all_kept <- do.call(computeSuperPopLDWindow, c(
    list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR",
         sample_keep = list(P2 = 12:1, P1 = 1:4, OTHER = 1L)),
    .covariance_fixture()
  ))
  expect_equal(all_kept$cormat, full$cormat, tolerance = 1e-12)
  expect_equal(all_kept$snplist, full$snplist)
})

test_that("sample_keep rejects missing populations and bad positions", {
  run <- function(keep) {
    do.call(computeSuperPopLDWindow, c(
      list(chr = 1, start_bp = 1, end_bp = 100, superpopulation_label = "EUR", sample_keep = keep),
      .covariance_fixture()
    ))
  }
  expect_error(run(list(P1 = 1:4)), "no entry for selected population P2")
  expect_error(run(list(P1 = 1:5, P2 = 1:12)), "out of range for population P1")
  expect_error(run(list(P1 = c(1L, 1L), P2 = 1:12)), "repeats a position")
  expect_error(run(list(P1 = integer(0), P2 = 1:12)), "keeps no individuals")
})
