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
