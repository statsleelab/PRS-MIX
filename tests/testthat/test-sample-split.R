.write_split_pop_desc <- function(path) {
  writeLines(c(
    "Population_Abbreviation\tNumber_of_Subject\tSuper_Population\tPopulation_Description",
    "P1\t10\tEUR\tone",
    "P2\t7\tEUR\ttwo",
    "P3\t9\tASN\tthree",
    "P4\t5\tAFR\tfour"
  ), path)
}

test_that("split_reference_samples splits every population disjointly and completely", {
  pop_desc <- tempfile(fileext = ".txt")
  .write_split_pop_desc(pop_desc)

  split <- split_reference_samples(pop_desc, c("EUR", "ASN"), seed = 1)

  expect_setequal(unique(split$population), c("P1", "P2", "P3"))
  sizes <- c(P1 = 10, P2 = 7, P3 = 9)
  for (p in names(sizes)) {
    sub <- split[split$population == p, ]
    # Every position appears exactly once.
    expect_identical(sort(sub$position), seq_len(sizes[[p]]))
  }
  # Largest-remainder sizes: 10 -> 5/3/2, 7 -> 3/2/2, 9 -> 5/2/2.
  tab <- table(split$population, factor(split$set, c("reference", "discovery", "target")))
  expect_equal(unname(tab["P1", ]), c(5, 3, 2))
  expect_equal(unname(tab["P2", ]), c(3, 2, 2))
  expect_equal(unname(tab["P3", ]), c(5, 2, 2))
})

test_that("split_reference_samples is reproducible by seed", {
  pop_desc <- tempfile(fileext = ".txt")
  .write_split_pop_desc(pop_desc)

  a <- split_reference_samples(pop_desc, "EUR", seed = 7)
  b <- split_reference_samples(pop_desc, "EUR", seed = 7)
  c <- split_reference_samples(pop_desc, "EUR", seed = 8)
  expect_identical(a, b)
  expect_false(identical(a$set, c$set))
})

test_that("split_reference_samples validates its inputs", {
  pop_desc <- tempfile(fileext = ".txt")
  .write_split_pop_desc(pop_desc)

  expect_error(split_reference_samples(pop_desc, "EUR", fractions = c(a = 0.5, b = 0.4), seed = 1),
               "summing to 1")
  expect_error(split_reference_samples(pop_desc, "XYZ", seed = 1), "No populations found")
  expect_error(split_reference_samples(pop_desc, "EUR"), "'seed'")
})

test_that("sample_keep_from_split returns sorted positions per population", {
  pop_desc <- tempfile(fileext = ".txt")
  .write_split_pop_desc(pop_desc)
  split <- split_reference_samples(pop_desc, c("EUR", "ASN"), seed = 3)

  keep <- sample_keep_from_split(split, "discovery")
  expect_setequal(names(keep), c("P1", "P2", "P3"))
  for (p in names(keep)) {
    expected <- sort(split$position[split$population == p & split$set == "discovery"])
    expect_identical(keep[[p]], expected)
  }
  expect_error(sample_keep_from_split(split, "nope"), "not found")
})
