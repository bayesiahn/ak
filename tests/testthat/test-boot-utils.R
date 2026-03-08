# Load the necessary package
library(testthat)

# Test that the function returns a vector of the correct length
test_that("generate_weights returns a vector of correct length", {
  N <- 10
  weights <- generate_weights(N)

  # Check if the length of the output matches N
  expect_equal(length(weights), N)
})

# Test that the weights sum to 1
test_that("generate_weights returns weights that sum to 1", {
  N <- 5
  weights <- generate_weights(N)

  # Check if the weights sum to 1
  expect_equal(sum(weights), 1)
})

# Test that all weights are non-negative
test_that("generate_weights returns non-negative weights", {
  N <- 7
  weights <- generate_weights(N)

  # Check that all weights are non-negative
  expect_true(all(weights >= 0))
})

# Test edge case with N = 1
test_that("generate_weights works for N = 1", {
  N <- 1
  weights <- generate_weights(N)

  # For N = 1, the single weight should be 1
  expect_equal(weights, 1)
})

# Test edge case with N = 0
test_that("generate_weights returns empty vector for N = 0", {
  N <- 0
  weights <- generate_weights(N)

  # Check that the output is an empty vector
  expect_equal(weights, numeric(0))
})


# Test for NULL unit_weights
test_that("check_unit_weights_validity defaults to equal weights when unit_weights is NULL", {
  posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
  unit_weights <- check_unit_weights_validity(posteriors, NULL)

  # Expected equal weights
  expected <- rep(1, nrow(posteriors))
  expect_equal(unit_weights, expected)
})

# Test for valid unit_weights
test_that("check_unit_weights_validity accepts valid unit_weights", {
  posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
  unit_weights <- c(1, 2)
  result <- check_unit_weights_validity(posteriors, unit_weights)

  # Should return the same unit_weights
  expect_equal(result, unit_weights)
})

# Test for mismatched unit_weights length
test_that("check_unit_weights_validity throws error for mismatched lengths", {
  posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
  unit_weights <- c(1, 2, 3)  # Incorrect length

  # Expect an error due to length mismatch
  expect_error(check_unit_weights_validity(posteriors, unit_weights),
               "Length of unit_weights must match the number of rows in posteriors.")
})

# Test for empty posteriors
test_that("check_unit_weights_validity works with empty posteriors", {
  posteriors <- matrix(numeric(0), nrow = 0, ncol = 2)
  unit_weights <- NULL

  # Expect an empty vector for unit_weights
  result <- check_unit_weights_validity(posteriors, unit_weights)
  expect_equal(result, numeric(0))
})

test_that("draw_multipliers returns correct length and values", {
  set.seed(123)
  N <- 100
  multipliers <- draw_multipliers(N)

  expect_equal(length(multipliers), N)  # Ensure correct length
  expect_true(all(multipliers %in% c(-1, 1)))  # Ensure only -1 and 1 are present
})

test_that("draw_multipliers returns a balanced distribution on large N", {
  set.seed(42)
  N <- 10000
  multipliers <- draw_multipliers(N)

  # Check approximate balance between -1 and 1
  expect_true(abs(sum(multipliers)) < N * 0.05)  # Should be close to zero
})

test_that("draw_multipliers handles edge cases correctly", {
  expect_equal(draw_multipliers(1) %in% c(-1, 1), TRUE)  # Single sample test
  expect_equal(length(draw_multipliers(0)), 0)  # Should return an empty vector
})

