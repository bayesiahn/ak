library(testthat)

# Unit tests
test_that("cbind_list_of_matrices works for list of matrices", {
  m1 <- matrix(1:4, nrow = 2)
  m2 <- matrix(5:8, nrow = 2)

  result <- cbind_list_of_matrices(list(m1, m2))
  expect_equal(result, cbind(m1, m2))
})

test_that("cbind_list_of_matrices works for list of numeric vectors", {
  v1 <- c(1, 2)
  v2 <- c(3, 4)

  result <- cbind_list_of_matrices(list(v1, v2))
  expect_equal(result, unname(cbind(v1, v2)))
})

test_that("cbind_list_of_matrices works for mixed matrices and vectors", {
  m1 <- matrix(1:4, nrow = 2)
  v2 <- c(5, 6)

  result <- cbind_list_of_matrices(list(m1, v2))
  expect_equal(result, unname(cbind(m1, v2)))
})

test_that("cbind_list_of_matrices works for nested list of matrices", {
  nested_list <- list(
    list(matrix(1:2, 2, 1), matrix(3:4, 2, 1)),
    list(matrix(5:6, 2, 1), matrix(7:8, 2, 1))
  )

  result <- cbind_list_of_matrices(nested_list)
  expected <- cbind(
    cbind(matrix(1:2, 2, 1), matrix(3:4, 2, 1)),
    cbind(matrix(5:6, 2, 1), matrix(7:8, 2, 1))
  )

  expect_equal(result, expected)
})

test_that("cbind_list_of_matrices throws error for invalid input", {
  invalid_list <- list("invalid", matrix(1:4, nrow=2))
  expect_error(cbind_list_of_matrices(invalid_list),
               "Input must be a list of matrices or a list of lists of matrices.")
})
