library(testthat)

test_that("get_y_t_by_Y_t_index", {
  # Test 1: Discrete case, exact matches
  test_that("Discrete case works with exact matches", {
    Y_t <- list(values = c(1, 2, 3), type = 'discrete')
    y_t <- c(1, 2, 3, 3)
    expect_equal(get_y_t_by_Y_t_index(y_t, Y_t), c(1, 2, 3, 3))
  })

  # Test 2: Discrete case, missing values
  test_that("Discrete case handles missing values", {
    Y_t <- list(values = c(1, 2, 3), type = 'discrete')
    y_t <- c(1, 4, 3, 0)
    expect_equal(get_y_t_by_Y_t_index(y_t, Y_t), c(1, NA, 3, NA))
  })

  # Test 3: Continuous case, exact and inexact matches
  test_that("Continuous case works with exact and inexact matches", {
    Y_t <- list(values = c(0.1, 0.3, 0.5), type = 'continuous')
    y_t <- c(0, 0.2, 0.4, 0.5)
    expect_equal(get_y_t_by_Y_t_index(y_t, Y_t), c(1, 2, 3, 3))
  })

  # Test 4: Continuous case, edge cases and out-of-bounds
  test_that("Continuous case handles edge cases and out-of-bounds values", {
    Y_t <- list(values = c(0.1, 0.3, 0.5), type = 'continuous')
    y_t <- c(0, 0.6, 0.1, 0.5)
    expect_equal(get_y_t_by_Y_t_index(y_t, Y_t), c(1, 4, 1, 3))
  })
})


test_that("get_transition_matrix", {
  # Test 1: Verify dimension of transition matrix
  test_that("Dimensions of transition matrix are correct", {
    y_matrix <- matrix(c(1, 2, 3, 2, 1, 3), ncol = 2)
    Y_t_now <- list(names = c(1, 2, 3), type = 'discrete')
    Y_t_forward <- list(names = c(1, 2, 3), type = 'discrete')
    result <- get_transition_matrix(y_matrix, Y_t_now, Y_t_forward)
    expect_equal(dim(result), c(3, 3), ignore_attr = TRUE)
  })

  # Test 2: Probabilities in rows sum to 1
  test_that("Row probabilities sum to 1", {
    y_matrix <- matrix(c(1, 2, 3, 2, 1, 3), ncol = 2)
    Y_t_now <- list(names = c(1, 2, 3), type = 'discrete')
    Y_t_forward <- list(names = c(1, 2, 3), type = 'discrete')
    result <- get_transition_matrix(y_matrix, Y_t_now, Y_t_forward)
    expect_true(all(rowSums(result) == 1))
  })

  # Test 3: Validate against a known transition matrix
  test_that("Correctly computes known transition matrix", {
    y_matrix <- matrix(c(1, 1, 1, 2, 2, 3), ncol = 2)
    Y_t_now <- list(names = c(1, 2, 3), type = 'discrete')
    Y_t_forward <- list(names = c(1, 2, 3), type = 'discrete')
    result <- get_transition_matrix(y_matrix, Y_t_now, Y_t_forward)
    known_result <- matrix(c(0, NaN, NaN, 2/3, NaN, NaN, 1/3, NaN, NaN), nrow = 3)
    expect_equal(result, known_result, ignore_attr = TRUE)
  })

  # Test 4: Handles cases with no transitions
  test_that("Handles no transitions correctly", {
    y_matrix <- matrix(c(1, 1, 1, 1, 1, 1), ncol = 2)
    Y_t_now <- list(names = c(1, 2, 3), type = 'discrete')
    Y_t_forward <- list(names = c(1, 2, 3), type = 'discrete')
    result <- get_transition_matrix(y_matrix, Y_t_now, Y_t_forward)
    expect_equal(result[1, 1], 1, ignore_attr = TRUE)
    expect_equal(sum(result[-1, ]), NaN, ignore_attr = TRUE)
  })

  # Test 5: Check if the function handles continuous type in Y_t
  test_that("Handles continuous type correctly", {
    y_matrix <- matrix(c(0.1, 0.3, 0.5, 0.3, 0.5, 0.1), ncol = 2)
    Y_t_now <- list(names = c(0.1, 0.3, 0.5), type = 'continuous')
    Y_t_forward <- list(names = c(0.1, 0.3, 0.5), type = 'continuous')
    result <- get_transition_matrix(y_matrix, Y_t_now, Y_t_forward)
    expect_true(all(rowSums(result, na.rm = TRUE) == 1))
  })

  # Test 6: Validate for non-square transition matrix
  test_that("Correctly generates non-square transition matrix", {
    y_matrix <- cbind(c(1, 2, 2), c(1, 3, 2))
    Y_t_now <- list(names = c(1, 2), type = 'discrete')
    Y_t_forward <- list(names = c(1, 2, 3), type = 'discrete')
    result <- get_transition_matrix(y_matrix, Y_t_now, Y_t_forward)
    expect_equal(dim(result), c(2, 3), ignore_attr = TRUE)
    expect_equal(result[1,], c(1, 0, 0), ignore_attr = TRUE)
    expect_equal(result[2,], c(0, 0.5, 0.5), ignore_attr = TRUE)
  })

  # Test 7: Validate handling of zero transitions for some states
  test_that("Handles zero transitions for some states correctly", {
    y_matrix <- matrix(c(1, 2, 1, 1, 2, 2), ncol = 2)
    Y_t_now <- list(names = c(1, 2, 3), type = 'discrete')
    Y_t_forward <- list(names = c(1, 2), type = 'discrete')
    result <- get_transition_matrix(y_matrix, Y_t_now, Y_t_forward)
    expect_true(is.nan(result[3, 1]))
    expect_true(is.nan(result[3, 2]))
  })

})


test_that("apply_transition_matrices computes correct PMF for a single matrix", {
  initial_pmf <- c(0.5, 0.3, 0.2)
  transition_matrices <- list(
    matrix(c(0.8, 0.1, 0.1,
             0.2, 0.6, 0.2,
             0.3, 0.3, 0.4), nrow = 3, byrow = TRUE)
  )

  result <- apply_transition_matrices(initial_pmf, transition_matrices)

  # Manually compute expected result
  expected <- initial_pmf %*% transition_matrices[[1]]
  expected <- as.numeric(expected / sum(expected))

  expect_equal(result, expected)
})

test_that("apply_transition_matrices computes correct PMF for multiple matrices", {
  initial_pmf <- c(0.5, 0.3, 0.2)
  transition_matrices <- list(
    matrix(c(0.8, 0.1, 0.1,
             0.2, 0.6, 0.2,
             0.3, 0.3, 0.4), nrow = 3, byrow = TRUE),
    matrix(c(0.7, 0.2, 0.1,
             0.2, 0.5, 0.3,
             0.4, 0.4, 0.2), nrow = 3, byrow = TRUE)
  )

  result <- apply_transition_matrices(initial_pmf, transition_matrices)

  # Manually compute expected result
  intermediate_pmf <- initial_pmf %*% transition_matrices[[1]]
  expected <- intermediate_pmf %*% transition_matrices[[2]]
  expected <- as.numeric(expected / sum(expected))

  expect_equal(result, expected)
})

test_that("apply_transition_matrices normalizes initial PMF if not normalized", {
  initial_pmf <- c(2, 3, 5)  # Not normalized
  transition_matrices <- list(
    matrix(c(0.8, 0.1, 0.1,
             0.2, 0.6, 0.2,
             0.3, 0.3, 0.4), nrow = 3, byrow = TRUE)
  )

  result <- apply_transition_matrices(initial_pmf, transition_matrices)

  # Normalize the initial PMF
  normalized_pmf <- initial_pmf / sum(initial_pmf)
  expected <- normalized_pmf %*% transition_matrices[[1]]
  expected <- as.numeric(expected / sum(expected))

  expect_equal(result, expected)
})

test_that("apply_transition_matrices handles edge case with empty transition_matrices", {
  initial_pmf <- c(0.5, 0.3, 0.2)
  transition_matrices <- list()  # No transition matrices

  result <- apply_transition_matrices(initial_pmf, transition_matrices)

  # Expected result should be the initial PMF, as no transitions are applied
  expect_equal(result, initial_pmf)
})

test_that("apply_transition_matrices throws error for invalid PMF", {
  initial_pmf <- c(-0.5, 0.3, 1.2)  # Invalid PMF (negative and sum > 1)
  transition_matrices <- list(
    matrix(c(0.8, 0.1, 0.1,
             0.2, 0.6, 0.2,
             0.3, 0.3, 0.4), nrow = 3, byrow = TRUE)
  )

  expect_error(
    apply_transition_matrices(initial_pmf, transition_matrices),
    "initial_pmf must be a numeric vector with non-negative values."
  )
})

test_that("apply_transition_matrices throws error for invalid transition_matrices", {
  initial_pmf <- c(0.5, 0.3, 0.2)
  transition_matrices <- list(
    c(0.8, 0.1, 0.1)  # Not a matrix
  )

  expect_error(
    apply_transition_matrices(initial_pmf, transition_matrices),
    "transition_matrices must be a list of matrices."
  )
})

test_that("apply_transition_matrices throws error for dimension mismatch", {
  initial_pmf <- c(0.5, 0.3, 0.2)
  transition_matrices <- list(
    matrix(c(0.8, 0.1,
             0.2, 0.6), nrow = 2, byrow = TRUE)  # Mismatch in dimensions
  )

  expect_error(
    apply_transition_matrices(initial_pmf, transition_matrices),
    "Each transition matrix must have a number of columns equal to the length of the PMF."
  )
})
