library(testthat)

test_that("bootstrap_transition_models_from_wide_matrix generates bootstrapped transition models", {
  y <- matrix(c(1, 2, 2,
                1, 1, 2), nrow = 2)
  g <- c(0, 3)
  opts <- list(bootstrap_counts = 10, mc_cores = 2)
  J <- 2

  result <- bootstrap_transition_models_from_wide_matrix(J, y, g, opts = opts)

  # Check the number of bootstrapped models
  expect_equal(length(result), opts$bootstrap_counts)

  # Check structure of one bootstrapped model
  expect_true(is.list(result[[1]]))
  expect_equal(length(result[[1]]$priors), J)
  expect_equal(length(result), opts$bootstrap_counts)
})

test_that("estimate_transition_model_from_wide_matrix generates bootstrapped transition models", {
  y <- matrix(c(1, 2, 2,
                1, 1, 2), nrow = 2)
  g <- c(0, 3)
  opts <- list(bootstrap_counts = 10, mc_cores = 2, se_method = "bootstrap")
  J <- 2

  result <- estimate_transition_model_from_wide_matrix(J, y, g, opts = opts)

  # Check the number of bootstrapped models
  expect_equal(length(result$bootstrap_estimates), opts$bootstrap_counts)

  # Check structure of one bootstrapped model
  expect_true(is.list(result$bootstrap_estimates[[1]]))
  expect_equal(length(result$bootstrap_estimates[[1]]$priors), J)
  expect_equal(length(result$bootstrap_estimates), opts$bootstrap_counts)
})
