# Load the necessary package
library(testthat)

# Test with equal weights (default case)
test_that("M_step_priors computes priors with equal weights", {
  posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
  result <- M_step_priors(posteriors)

  # Expected priors
  expected <- c((0.2 + 0.3) / 1.0, (0.4 + 0.1) / 1.0)
  expect_equal(result, expected)
})

# Test with provided unit weights
test_that("M_step_priors computes priors with provided unit weights", {
  posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
  unit_weights <- c(1, 2)
  result <- M_step_priors(posteriors, unit_weights)

  # Expected weighted priors
  expected <- c((0.2 + 0.3 * 2), (0.4 + 0.1 * 2))
  expected <- expected / sum(expected)
  expect_equal(result, expected)
})

test_that("M_step_priors computes priors with provided unit weights, J =1 ", {
  posteriors <- matrix(c(0.2, 0.3), ncol = 1)
  unit_weights <- c(1, 2)
  result <- M_step_priors(posteriors, unit_weights)

  # Expected weighted priors
  expected <- 1
  expect_equal(result, expected)
})

# Test with unequal weights
test_that("M_step_priors works with arbitrary unit weights", {
  posteriors <- matrix(c(0.1, 0.4, 0.2, 0.3, 0.3, 0.2), ncol = 3)
  unit_weights <- c(2, 1)
  result <- M_step_priors(posteriors, unit_weights)

  # Expected weighted priors
  weighted_posteriors <- posteriors * unit_weights
  expected <- colSums(weighted_posteriors) / sum(colSums(weighted_posteriors))
  expect_equal(result, expected)
})

# Test with no unit_weights (equal weights implied)
test_that("M_step_priors assumes equal weights when unit_weights is NULL", {
  posteriors <- matrix(c(0.1, 0.4, 0.2, 0.3, 0.3, 0.2), ncol = 3)
  result <- M_step_priors(posteriors)

  # Expected priors for equal weights
  expected <- colSums(posteriors) / sum(colSums(posteriors))
  expect_equal(result, expected)
})

# Test edge case with unit_weights of all zeros
test_that("M_step_priors handles unit_weights with all zeros", {
  posteriors <- matrix(c(0.1, 0.4, 0.2, 0.3), ncol = 2)
  unit_weights <- c(0, 0)

  expect_equal(M_step_priors(posteriors, unit_weights), c(NaN, NaN))
})

# Test with mismatched dimensions
test_that("M_step_priors throws error for mismatched unit_weights and posteriors", {
  posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
  unit_weights <- c(1, 2, 3)  # Incorrect length

  # Expect an error due to dimension mismatch
  expect_error(M_step_priors(posteriors, unit_weights),
               "Length of unit_weights must match the number of rows in posteriors.")
})

# Test for equal weights (default case)
test_that("M_step_initial_outcome_dists computes distributions with equal weights", {
  posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
  indicators <- rbind(c(1, 0), c(0, 1))

  result <- M_step_initial_outcome_dists(posteriors, indicators)

  # Expected distributions for each latent group
  expected <- list(
    c(0.2, 0.3)/0.5,
    c(0.4, 0.1)/0.5
  )
  expect_equal(result, expected)
})

# Test for provided unit weights
test_that("M_step_initial_outcome_dists computes distributions with provided unit weights", {
  posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
  indicators <- rbind(c(1, 0), c(0, 1))
  unit_weights <- c(1, 2)

  result <- M_step_initial_outcome_dists(posteriors, indicators, unit_weights)

  # Expected distributions with unit weights applied
  expected <- list(
    c(0.2 * 1 + 0 * 2, 0.3 * 2 + 0 * 1)/sum(c(0.2 * 1 + 0 * 2, 0.3 * 2 + 0 * 1)),
    c(0.4 * 1 + 0 * 2, 0.1 * 2 + 0 * 1)/sum(c(0.4 * 1 + 0 * 2, 0.1 * 2 + 0 * 1))
  )
  expect_equal(result, expected)
})

# Test with single latent group
test_that("M_step_initial_outcome_dists works with single latent group", {
  posteriors <- matrix(c(0.2, 0.3), ncol = 1)
  indicators <- rbind(c(1, 0), c(0, 1))
  unit_weights <- c(1, 2)

  result <- M_step_initial_outcome_dists(posteriors, indicators, unit_weights)

  # Expected distribution for single group with weights
  expected <- list(
    c(0.2 * 1 + 0 * 2, 0.3 * 2 + 0 * 1) / sum(c(0.2 * 1 + 0 * 2, 0.3 * 2 + 0 * 1))
  )
  expect_equal(result, expected)
})

# Test with no unit_weights (equal weights implied)
test_that("M_step_initial_outcome_dists assumes equal weights when unit_weights is NULL", {
  posteriors <- matrix(c(0.1, 0.4, 0.2, 0.3, 0.3, 0.2), ncol = 3)
  indicators <- rbind(c(1, 0), c(0, 1))

  result <- M_step_initial_outcome_dists(posteriors, indicators)

  # Expected distributions assuming equal weights
  expected <- list(
    c(0.1, 0.4)/0.5,
    c(0.2, 0.3)/0.5,
    c(0.3, 0.2)/0.5
  )
  expect_equal(result, expected)
})

# Test edge case with empty indicators
test_that("M_step_initial_outcome_dists returns NULL for empty indicators", {
  posteriors <- matrix(c(0.1, 0.4), ncol = 1)
  indicators <- rbind()

  result <- M_step_initial_outcome_dists(posteriors, indicators)

  expect_null(result)
})

# Test for mismatched dimensions
test_that("M_step_initial_outcome_dists throws error for mismatched unit_weights and posteriors", {
  posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
  indicators <- rbind(c(1, 0), c(0, 1))
  unit_weights <- c(1, 2, 3)  # Incorrect length

  expect_error(
    M_step_initial_outcome_dists(posteriors, indicators, unit_weights),
    "Length of unit_weights must match the number of rows in posteriors."
  )
})


test_that("M_step_Ps computes transition matrices with equal weights, J = 1", {
  posteriors <- matrix(c(1, 1, 1), ncol = 1)
  indicators <- list(
    list(matrix(c(1, 0, 0, 0), 2, 2), matrix(c(0, 0, 1, 0), 2, 2), matrix(c(0, 1, 0, 0), 2, 2)),
    list(matrix(c(0, 1, 0, 0), 2, 2), matrix(c(0, 0, 0, 1), 2, 2), matrix(c(1, 0, 0, 0), 2, 2))
  )

  result <- M_step_Ps(posteriors, indicators)

  # Expected transition matrices for each group
  expected <- list(list(
    matrix(c(0.5, 1, 0.5, 0), 2, 2),
    matrix(c(1, 0.5, 0, 0.5), 2, 2)
  ))
  expect_equal(result, expected)
})


test_that("M_step_Ps computes transition matrices with unequal weights, J = 1", {
  posteriors <- matrix(c(3, 1, 1), ncol = 1)
  indicators <- list(
    list(matrix(c(1, 0, 0, 0), 2, 2), matrix(c(0, 0, 1, 0), 2, 2), matrix(c(0, 1, 0, 0), 2, 2)),
    list(matrix(c(0, 1, 0, 0), 2, 2), matrix(c(0, 0, 0, 1), 2, 2), matrix(c(1, 0, 0, 0), 2, 2))
  )

  result <- M_step_Ps(posteriors, indicators)

  # Expected transition matrices for each group
  expected <- list(list(
    matrix(c(0.75, 1, 0.25, 0), 2, 2),
    matrix(c(1, 0.75, 0, 0.25), 2, 2)
  ))
  expect_equal(result, expected)
})

test_that("M_step_Ps computes transition matrices with equal weights, J > 1", {
  posteriors <- matrix(c(1, 1, 1, 1, 1, 1), ncol = 2)
  indicators <- list(
    list(matrix(c(1, 0, 0, 0), 2, 2), matrix(c(0, 0, 1, 0), 2, 2), matrix(c(0, 1, 0, 0), 2, 2)),
    list(matrix(c(0, 1, 0, 0), 2, 2), matrix(c(0, 0, 0, 1), 2, 2), matrix(c(1, 0, 0, 0), 2, 2))
  )

  result <- M_step_Ps(posteriors, indicators)

  # Expected transition matrices for each group
  expected <- list(
    list(
      matrix(c(0.5, 1, 0.5, 0), 2, 2),
      matrix(c(1, 0.5, 0, 0.5), 2, 2)
    ),
    list(
      matrix(c(0.5, 1, 0.5, 0), 2, 2),
      matrix(c(1, 0.5, 0, 0.5), 2, 2)
    )
  )
  expect_equal(result, expected)
})

test_that("M_step_Ps computes transition matrices with unequal weights, J > 1", {
  posteriors <- matrix(c(3, 1, 1,
                         3, 2, 1), ncol = 2)
  indicators <- list(
    list(matrix(c(1, 0, 0, 0), 2, 2), matrix(c(0, 0, 1, 0), 2, 2), matrix(c(0, 1, 0, 0), 2, 2)),
    list(matrix(c(0, 1, 0, 0), 2, 2), matrix(c(0, 0, 0, 1), 2, 2), matrix(c(1, 0, 0, 0), 2, 2))
  )

  result <- M_step_Ps(posteriors, indicators)

  # Expected transition matrices for each group
  expected <- list(
    list(
      matrix(c(0.75, 1, 0.25, 0), 2, 2),
      matrix(c(1, 0.75, 0, 0.25), 2, 2)
    ),
    list(
      matrix(c(0.6, 1, 0.4, 0), 2, 2),
      matrix(c(1, 0.6, 0, 0.4), 2, 2)
    )
  )
  expect_equal(result, expected)
})

