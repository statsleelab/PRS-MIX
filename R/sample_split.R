#' Split reference individuals into disjoint sets, within each population
#'
#' Assigns every individual of the selected superpopulations to one of several
#' disjoint sets (by default 50% reference / 25% discovery / 25% target),
#' separately within each population so every set keeps the population mix of
#' the full superpopulation. Individuals are identified by their 1-based
#' position within their population's genotype string, the only identifier the
#' 33KG reference files carry.
#'
#' Set sizes within a population are the floor of `fractions * n`, with the
#' remaining individuals going to the sets with the largest remainders (ties
#' broken by the order of `fractions`), so sizes are deterministic and sum to n.
#' Which individuals go where is random, controlled by `seed`.
#'
#' Relatedness cannot be checked from the reference files, so relatives may
#' fall into different sets.
#'
#' @param reference_pop_desc_file Reference panel population description file.
#' @param superpopulation_labels Superpopulations to split (e.g. `c("EUR", "ASN")`).
#' @param fractions Named numeric vector of set fractions summing to 1.
#' @param seed Integer RNG seed.
#' @return A data frame with one row per individual: `population`,
#'   `super_population`, `position` (1-based within the population) and `set`.
#' @export
split_reference_samples <- function(reference_pop_desc_file,
                                    superpopulation_labels,
                                    fractions = c(reference = 0.5, discovery = 0.25, target = 0.25),
                                    seed) {
  if (!is.numeric(fractions) || is.null(names(fractions)) || any(!nzchar(names(fractions))) ||
      anyDuplicated(names(fractions)) || any(!is.finite(fractions)) || any(fractions <= 0) ||
      abs(sum(fractions) - 1) > 1e-8) {
    stop("'fractions' must be a named vector of positive fractions summing to 1.", call. = FALSE)
  }
  if (missing(seed) || length(seed) != 1 || !is.finite(seed)) {
    stop("'seed' must be a single number.", call. = FALSE)
  }

  pop_desc <- .read_pop_desc(reference_pop_desc_file)
  size_col <- grep("^Number_of_Subject", names(pop_desc), value = TRUE)
  if (length(size_col) != 1) {
    stop("reference_pop_desc_file must have one Number_of_Subject(s) column.", call. = FALSE)
  }
  labels <- toupper(superpopulation_labels)
  pops <- pop_desc[toupper(pop_desc$Super_Population) %in% labels, , drop = FALSE]
  missing_labels <- setdiff(labels, toupper(pops$Super_Population))
  if (length(missing_labels) > 0) {
    stop(sprintf("No populations found for: %s", paste(missing_labels, collapse = ", ")), call. = FALSE)
  }

  set.seed(seed)
  parts <- lapply(seq_len(nrow(pops)), function(i) {
    n <- as.integer(pops[[size_col]][i])
    sizes <- .split_sizes(n, fractions)
    data.frame(
      population = pops$Population_Abbreviation[i],
      super_population = toupper(pops$Super_Population[i]),
      position = sample.int(n),
      set = rep(names(fractions), sizes),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, parts)
  out <- out[order(out$super_population, out$population, out$position), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Largest-remainder allocation of n individuals to fractions.
.split_sizes <- function(n, fractions) {
  raw <- fractions * n
  sizes <- floor(raw)
  extra <- n - sum(sizes)
  if (extra > 0) {
    rem <- raw - sizes
    give <- order(-rem, seq_along(rem))[seq_len(extra)]
    sizes[give] <- sizes[give] + 1
  }
  as.integer(sizes)
}

#' Turn one set of a sample split into a `sample_keep` list
#'
#' @param split Output of \code{\link{split_reference_samples}}.
#' @param set Name of the set to keep (e.g. `"reference"`).
#' @return Named list (population -> sorted 1-based positions), as accepted by
#'   the `sample_keep` argument of \code{\link{computeSuperPopLDWindow}},
#'   \code{\link{readSuperPopGenotypes}} and
#'   \code{\link{generate_prscs_ld_panels}}.
#' @export
sample_keep_from_split <- function(split, set) {
  if (!all(c("population", "position", "set") %in% names(split))) {
    stop("'split' must have columns population, position and set.", call. = FALSE)
  }
  if (!set %in% split$set) {
    stop(sprintf("Set '%s' not found in split.", set), call. = FALSE)
  }
  sub <- split[split$set == set, , drop = FALSE]
  lapply(split(sub$position, sub$population), function(p) sort(as.integer(p)))
}
