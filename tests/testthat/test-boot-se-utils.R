library(testthat)

# Test scalar inputs
test_that("elementwise_sd handles scalar inputs", {
  values <- list(1, 2, 3)
  result <- elementwise_sd(values)
  expected <- sqrt(var(c(1, 2, 3)))
  expect_equal(result, expected)
})

# Test vector inputs
test_that("elementwise_sd handles vector inputs", {
  values <- list(c(1, 2), c(2, 3), c(3, 4))
  result <- elementwise_sd(values)

  # Expected standard errors for each element
  expected <- c(sqrt(var(c(1, 2, 3))), sqrt(var(c(2, 3, 4))))
  expect_equal(result, expected)
})

# Test matrix inputs
test_that("elementwise_sd handles matrix inputs", {
  values <- list(matrix(1:4, 2, 2), matrix(2:5, 2, 2))
  result <- elementwise_sd(values)

  # Expected element-wise standard errors
  expected <- matrix(apply(array(c(1, 2, 3, 4, 2, 3, 4, 5), dim = c(2, 2, 2)), c(1, 2), sd), 2, 2)
  expect_equal(result, expected)
})

# Test nested list of matrices
test_that("elementwise_sd handles nested list of matrices", {
  values <- list(
    list(matrix(1, 2, 2), matrix(1, 2, 2)),
    list(matrix(2, 2, 2), matrix(2, 2, 2))
  )
  result <- elementwise_sd(values)

  # Expected element-wise standard errors for each nested matrix
  expected <- list(
    matrix(c(sqrt(0.5), sqrt(0.5), sqrt(0.5), sqrt(0.5)), 2, 2),
    matrix(c(sqrt(0.5), sqrt(0.5), sqrt(0.5), sqrt(0.5)), 2, 2)
  )
  expect_equal(result, expected)
})

# Test empty input
test_that("elementwise_sd handles empty input", {
  values <- list()
  result <- elementwise_sd(values)
  expect_equal(result, list())
})

test_that("elementwise_quantile computes quantiles for scalar inputs", {
  values <- list(1, 2, 3)
  result <- elementwise_quantile(values, probs = c(0.25, 0.5, 0.75))

  # Expected quantiles for scalar values
  expected <- as.matrix(quantile(c(1, 2, 3), probs = c(0.25, 0.5, 0.75)))
  expect_equal(result, expected)
})

test_that("elementwise_quantile computes quantiles for vector inputs", {
  values <- list(c(1, 2), c(2, 3), c(3, 4))
  result <- elementwise_quantile(values, probs = c(0.25, 0.5, 0.75))

  # Expected quantiles for each element
  expected <- apply(do.call(rbind, values), 2, quantile, probs = c(0.25, 0.5, 0.75))
  expect_equal(result, expected)
})

test_that("elementwise_quantile computes quantiles for matrix inputs", {
  values <- list(matrix(1:4, 2, 2), matrix(2:5, 2, 2))
  result <- elementwise_quantile(values, probs = c(0.25, 0.5, 0.75))

  # Expected element-wise quantiles
  array_values <- array(c(1, 2, 3, 4, 2, 3, 4, 5), dim = c(2, 2, 2))
  expected <- apply(array_values, c(1, 2), quantile, probs = c(0.25, 0.5, 0.75))
  expect_equal(result, expected)
})

test_that("elementwise_quantile computes quantiles for nested lists of matrices", {
  values <- list(
    list(matrix(1:4, 2, 2), matrix(5:8, 2, 2)),
    list(matrix(2:5, 2, 2), matrix(6:9, 2, 2))
  )
  result <- elementwise_quantile(values, probs = c(0.25, 0.5, 0.75))

  # Expected element-wise quantiles for nested structures
  expected <- list(
      array(c(1.25, 1.50, 1.75, 2.25, 2.50, 2.75, 3.25, 3.50, 3.75, 4.25, 4.50, 4.75),
            dim = c(3, 2, 2), dimnames = list(c("25%", "50%", "75%"), NULL, NULL)),
      array(c(5.25, 5.50, 5.75, 6.25, 6.50, 6.75, 7.25, 7.50, 7.75, 8.25, 8.50, 8.75),
            dim = c(3, 2, 2), dimnames = list(c("25%", "50%", "75%"), NULL, NULL)
    )
  )

  expect_equal(result, expected)
})

test_that("elementwise_quantile handles empty input", {
  values <- list()
  result <- elementwise_quantile(values, probs = c(0.25, 0.5, 0.75))
  expect_equal(result, list())
})

test_that("elementwise_quantile handles single matrix input", {
  values <- list(matrix(1:4, 2, 2))
  result <- elementwise_quantile(values, probs = c(0.25, 0.5, 0.75))

  # Expected element-wise quantiles
  expected <- apply(array(c(1, 2, 3, 4), dim = c(2, 2, 1)), c(1, 2), quantile, probs = c(0.25, 0.5, 0.75))
  expect_equal(result, expected)
})
