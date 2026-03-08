# Unit Tests for generate_Y_t function
test_that("generate_Y_t works for discrete variables", {
  y_t_discrete <- c(1, 2, 3, 3)
  result <- generate_Y_t(y_t_discrete)
  expect_equal(result$type, "discrete")
  expect_equal(result$name, (sort(unique(y_t_discrete))))
  expect_equal(result$values, c(0, 0, 1))
})

test_that("get_interval_names works as expected", {
  expect_equal(get_interval_names(c(0, 0.1, 0.3, 0.5)), c("[0.00, 0.10)", "[0.10, 0.30)", "[0.30, 0.50]"))
  expect_equal(get_interval_names(c(0, 1)), c("[0.00, 1.00]"))
  expect_equal(get_interval_names(numeric(0)), character(0))
  expect_equal(get_interval_names(c(1, 1)), "1.00")
})

test_that("get_conditional_state_matrix works with lags = 1", {
  y <- matrix(1:6, nrow = 2, byrow = TRUE)
  result <- get_conditional_state_matrix(y, lags = 1)
  expected <- matrix(c(NA, "1_2", "2_3", NA, "4_5", "5_6"), nrow = 2, byrow = TRUE)
  expect_equal(result, expected)
})

test_that("get_conditional_state_matrix works with lags = 2", {
  y <- matrix(1:6, nrow = 2, byrow = TRUE)
  result <- get_conditional_state_matrix(y, lags = 2)
  expected <- matrix(c(NA, NA, "1_2_3", NA, NA, "4_5_6"), nrow = 2, byrow = TRUE)
  expect_equal(result, expected)
})

test_that("get_conditional_state_matrix works with lags = 0", {
  y <- matrix(1:6, nrow = 2, byrow = TRUE)
  result <- get_conditional_state_matrix(y, lags = 0)
  expected <- y
  expect_equal(result, expected)
})

test_that("get_conditional_state_matrix returns NA for insufficient lags", {
  y <- matrix(1:6, nrow = 2, byrow = TRUE)
  result <- get_conditional_state_matrix(y, lags = 3)
  expected <- matrix(NA, nrow = 2, ncol = 3)
  expect_equal(result, expected)
})

test_that("get_conditional_state_matrix handles edge case of single row matrix", {
  y <- matrix(1:3, nrow = 1)
  result <- get_conditional_state_matrix(y, lags = 1)
  expected <- matrix(c(NA, "1_2", "2_3"), nrow = 1)
  expect_equal(result, expected)
})

test_that("get_conditional_state_matrix handles edge case of single column matrix", {
  y <- matrix(1:3, ncol = 1)
  result <- get_conditional_state_matrix(y, lags = 1)
  expected <- matrix(NA, ncol = 1, nrow = 3)
  expect_equal(result, expected)
})
