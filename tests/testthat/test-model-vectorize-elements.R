library(testthat)

test_that("vectorize_pmf correctly vectorizes PMF", {
  expect_equal(vectorize_pmf(c(0.3, 0.7)), c(0.3))
  expect_equal(vectorize_pmf(c(1.0)), numeric(0))
  expect_equal(vectorize_pmf(c(0.4, 0.6), minimal = FALSE), c(0.4, 0.6))
})

test_that("vectorize_pmfs applies vectorization to a list of PMFs", {
  pmfs <- list(c(0.3, 0.7), c(0.2, 0.8), c(0.5))
  result <- vectorize_pmfs(pmfs)
  expect_equal(result, list(c(0.3), c(0.2), numeric(0)))
})

test_that("vectorize_transition_matrix correctly vectorizes transition matrices", {
  matrix <- matrix(c(0.6, 0.4), nrow = 1)
  expect_equal(vectorize_transition_matrix(matrix), c(0.6))

  matrix2 <- matrix(c(0.6, 0.4, 0.2, 0.8), nrow = 2)
  expect_equal(vectorize_transition_matrix(matrix2), c(0.6, 0.4))
})

test_that("vectorize_transition_matrices applies vectorization to a list of matrices", {
  matrices <- list(matrix(c(0.6, 0.4), nrow = 1),
                   matrix(c(0.7, 0.3, 0.5, 0.5), nrow = 2))
  result <- vectorize_transition_matrices(matrices)
  expect_equal(result, list(c(0.6), c(0.7, 0.3)))
})

test_that("vectorization functions handle minimal = FALSE correctly", {
  expect_equal(vectorize_pmf(c(0.3, 0.7), minimal = FALSE), c(0.3, 0.7))

  matrix <- matrix(c(0.6, 0.4), nrow = 1)
  expect_equal(vectorize_transition_matrix(matrix, minimal = FALSE), c(0.6, 0.4))

  pmfs <- c(c(0.3, 0.7), c(0.2, 0.8))
  expect_equal(vectorize_pmfs(pmfs, minimal = FALSE), pmfs)
})
