library(testthat)

test_that("unvectorize_pmf correctly reconstructs PMF", {
  expect_equal(unvectorize_pmf(c(0.3, 0.5)), c(0.3, 0.5, 0.2))
  expect_equal(unvectorize_pmf(numeric(0)), c(1.0))
  expect_equal(unvectorize_pmf(c(0.6, 0.4), minimal = FALSE), c(0.6, 0.4))
})
test_that("unvectorize_transition_matrix correctly reconstructs transition matrix", {
  vectorized <- c(0.7, 0.2, 0.1, 0.3)
  expected_matrix <- matrix(c(0.7, 0.2, 0.1, 0.1, 0.3, 0.6), nrow = 2, byrow = TRUE)

  result <- unvectorize_transition_matrix(vectorized, nrow = 2)
  expect_equal(result, expected_matrix)
})

test_that("unvectorize_transition_matrix handles minimal = FALSE", {
  vectorized <- c(0.6, 0.4, 0.5, 0.5)
  expected_matrix <- matrix(c(0.6, 0.4, 0.5, 0.5), nrow = 2, byrow = TRUE)

  result <- unvectorize_transition_matrix(vectorized, nrow = 2, minimal = FALSE)
  expect_equal(result, expected_matrix)
})

test_that("unvectorize_transition_matrix handles single-row case", {
  vectorized <- c(0.8, 0.2)
  expected_matrix <- matrix(c(0.8, 0.2), nrow = 1, byrow = TRUE)

  result <- unvectorize_transition_matrix(vectorized, nrow = 1, minimal = FALSE)
  expect_equal(result, expected_matrix)
})

test_that("vectorizing then unvectorizing PMF returns original PMF", {
  for (minimal in c(TRUE, FALSE)) {
    pmf <- c(0.3, 0.5, 0.2)
    vectorized <- vectorize_pmf(pmf, minimal = minimal)
    restored <- unvectorize_pmf(vectorized, minimal = minimal)

    expect_equal(restored, pmf)
  }
})

test_that("vectorizing then unvectorizing a single-value PMF returns the same", {
  for (minimal in c(TRUE, FALSE)) {
    pmf <- c(1.0)
    vectorized <- vectorize_pmf(pmf, minimal = minimal)
    restored <- unvectorize_pmf(vectorized, minimal = minimal)

    expect_equal(restored, pmf)
  }
})

test_that("vectorizing then unvectorizing transition matrix returns original", {
  for (nrow in 2:4) {
    for (minimal in c(TRUE, FALSE)) {
    transition_matrix <- generate_transition_model(outcomes = 1:nrow)$Ps_control[[1]][[1]]
    vectorized <- vectorize_transition_matrix(transition_matrix, minimal)
    restored <- unvectorize_transition_matrix(vectorized, nrow = nrow, minimal)

    expect_equal(restored, unname(transition_matrix))
    }
  }
})

test_that("vectorizing then unvectorizing a 1-row transition matrix returns original", {
  for (minimal in c(TRUE, FALSE)) {
    transition_matrix <- matrix(c(0.8, 0.2), nrow = 1, byrow = TRUE)
    vectorized <- vectorize_transition_matrix(transition_matrix, minimal)
    restored <- unvectorize_transition_matrix(vectorized, nrow = 1, minimal)

    expect_equal(restored, transition_matrix)
  }
})

test_that("unvectorize_transition_matrices correctly reconstructs multiple matrices", {
  vectorized <- c(0.7, 0.3, 0.4, 0.6, 0.1, 0.4, 0.5, 0.4, 0.5, 0.1)
  vectorized_minimal <- c(0.7, 0.4, 0.1, 0.4, 0.4, 0.5)
  nrows <- c(2, 2)  # Two 2x2 matrices

  expected_matrices <- list(
    matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE),
    matrix(c(0.1, 0.4, 0.5, 0.4, 0.5, 0.1), nrow = 2, byrow = TRUE)
  )

  result <- unvectorize_transition_matrices(vectorized, nrows, minimal = FALSE)
  result_minimal <- unvectorize_transition_matrices(vectorized_minimal, nrows, minimal = TRUE)
  expect_equal(result, expected_matrices)
  expect_equal(result_minimal, expected_matrices)
})

test_that("unvectorize_transition_matrices correctly reconstructs a single matrix", {
  vectorized <- c(0.8, 0.6)
  nrows <- c(2)  # One 2x2 matrix

  expected_matrix <- list(matrix(c(0.8, 0.2, 0.6, 0.4), nrow = 2, byrow = TRUE))

  result <- unvectorize_transition_matrices(vectorized, nrows, minimal = TRUE)
  expect_equal(result, expected_matrix)
})

test_that("unvectorize_transition_matrices handles minimal = FALSE correctly", {
  vectorized <- c(0.6, 0.4, 0.5, 0.5)
  nrows <- c(2)

  expected_matrices <- list(matrix(c(0.6, 0.4, 0.5, 0.5), nrow = 2, byrow = TRUE))

  result <- unvectorize_transition_matrices(vectorized, nrows, minimal = FALSE)
  expect_equal(result, expected_matrices)
})

test_that("unvectorize_transition_matrices handles a matrix with a single row", {
  vectorized <- c(0.8, 0.15)
  nrows <- c(1)  # One 1x3 matrix

  expected_matrix <- list(matrix(c(0.8, 0.15, 0.05), nrow = 1, byrow = TRUE))

  result <- unvectorize_transition_matrices(vectorized, nrows, minimal = TRUE)
  expect_equal(result, expected_matrix)
})

test_that("unvectorize_transition_matrices handles empty input", {
  expect_equal(unvectorize_transition_matrices(numeric(0), numeric(0)), list())
})

test_that("unvectorize_transition_matrices handles inconsistent row sizes", {
  vectorized <- c(0.6, 0.2, 0.1, 0.4, 0.4, 0.8, 0.15)
  nrows <- c(2, 3)  # One 2x3 and one 3x2 matrix

  expected_matrices <- list(
    matrix(c(0.6, 0.2, 0.2,
             0.1, 0.4, 0.5), nrow = 2, byrow = TRUE),
    matrix(c(0.4, 0.6,
             0.8, 0.2,
             0.15, 0.85), nrow = 3, byrow = TRUE)
  )

  result <- unvectorize_transition_matrices(vectorized, nrows, minimal = TRUE)
  expect_equal(result, expected_matrices)
})

test_that("unvectorize_pmfs correctly reconstructs multiple lists of PMFs", {
  vectorized_pmfs <- c(0.3, 0.5, 0.4, 0.4)
  list_length <- 2

  expected_pmfs <- list(
    c(0.3, 0.5, 0.2),
    c(0.4, 0.4, 0.2)
  )

  result <- unvectorize_pmfs(vectorized_pmfs, list_length)
  expect_equal(result, expected_pmfs)
})

test_that("unvectorize_pmfs correctly reconstructs a single list", {
  vectorized_pmfs <- c(0.6, 0.4)
  list_length <- 1

  expected_pmfs <- list(
    (c(0.6, 0.4, 0.0))
  )

  result <- unvectorize_pmfs(vectorized_pmfs, list_length)
  expect_equal(result, expected_pmfs)
})

test_that("unvectorize_pmfs handles minimal = FALSE", {
  vectorized_pmfs <- c(0.6, 0.4, 0.5, 0.5)
  list_length <- 2

  expected_pmfs <- list(
    (c(0.6, 0.4)),
    (c(0.5, 0.5))
  )

  result <- unvectorize_pmfs(vectorized_pmfs, list_length, minimal = FALSE)
  expect_equal(result, expected_pmfs)
})

test_that("unvectorize_pmfs handles empty input", {
  expect_equal(unvectorize_pmfs(numeric(0), list_length = 2), list(list(), list()))
})

test_that("unvectorize_pmfs handles inconsistent lengths", {
  vectorized_pmfs <- c(0.1, 0.5, 0.2, 0.4, 0.6, 0.35, 0.0)
  list_length <- 3

  expect_error(unvectorize_pmfs(vectorized_pmfs, list_length),
               "Length of vectorized_pmfs must be divisible by list_length")
})

test_that("unvectorize_transition_matrices_list correctly reconstructs multiple lists", {
  vectorized <- c(0.7, 0.2,
                  0.1, 0.3,
                  0.1, 0.7,
                  0.6, 0.4)
  nrows <- c(2, 2)
  list_length <- 2

  expected_matrices <- list(
    list(matrix(c(0.7, 0.2, 0.3, 0.8), ncol = 2),
              matrix(c(0.1, 0.3, 0.9, 0.7), ncol = 2)),
    list(matrix(c(0.1, 0.7, 0.9, 0.3), ncol = 2),
              matrix(c(0.6, 0.4, 0.4, 0.6), ncol = 2))
  )

  result <- unvectorize_transition_matrices_list(vectorized, nrows, list_length,
                                                 minimal = TRUE)
  expect_equal(result, expected_matrices)
})

test_that("unvectorize_transition_matrices_list correctly reconstructs a single list", {
  vectorized <- c(0.8,
                  0.6)
  nrows <- c(2)
  list_length <- 1

  expected_matrices <- list(
    list((matrix(c(0.8, 0.2, 0.6, 0.4), nrow = 2, byrow = TRUE)))
  )

  result <- unvectorize_transition_matrices_list(vectorized, nrows, list_length)
  expect_equal(result, expected_matrices)
})

test_that("unvectorize_transition_matrices_list handles minimal = FALSE", {
  vectorized <- c(0.6, 0.4, 0.5, 0.5)
  nrows <- c(2)
  list_length <- 1

  expected_matrices <- list(
    list(matrix(c(0.6, 0.4, 0.5, 0.5), nrow = 2, byrow = TRUE))
  )

  result <- unvectorize_transition_matrices_list(vectorized, nrows, list_length, minimal = FALSE)
  expect_equal(result, expected_matrices)
})

test_that("unvectorize_transition_matrices_list handles a matrix with a single row", {
  vectorized <- c(0.7, 0.1)
  nrows <- c(1)
  list_length <- 1

  expected_matrices <- list(
    list(matrix(c(0.7, 0.1, 0.2), nrow = 1, byrow = TRUE))
  )

  result <- unvectorize_transition_matrices_list(vectorized, nrows, list_length)
  expect_equal(result, expected_matrices)
})

test_that("unvectorize_transition_matrices_list handles empty input", {
  expect_equal(unvectorize_transition_matrices_list(numeric(0), nrows = numeric(0), list_length = 2), list(list(), list()))
})

test_that("unvectorize_transition_matrices_list handles inconsistent list_length", {
  vectorized <- c(0.6, 0.4, 0.3, 0.7, 0.8, 0.2, 0.1)
  nrows <- c(2, 3)  # One 2x2 and one 3x2 matrix
  list_length <- 2

  expect_error(unvectorize_transition_matrices_list(vectorized, nrows, list_length),
               "Length of vectorized_matrices_list must be divisible by list_length")
})
