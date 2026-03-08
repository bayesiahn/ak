set.seed(123456)

# Unit Tests for get_pmf function
test_that("get_pmf works for discrete variables", {
  y_t_discrete <- c(1, 1, 2, 2, 2, 3)
  Y_t_discrete <- list(values = c(1, 2, 3), type = 'discrete')
  result <- get_pmf(y_t_discrete, Y_t_discrete)
  expect_equal(result, c(2/6, 3/6, 1/6), ignore_attr = TRUE)
  expect_equal(names(result), get_Y_t_names(Y_t_discrete))
})

test_that("get_pmf works for continuous variables with grid_length", {
  y_t_continuous <- c(0.1, 0.2, 0.1, 0.3)
  Y_t_continuous <- list(values = c(0.1, 0.2, 0.3), type = 'continuous')
  result <- get_pmf(y_t_continuous, Y_t_continuous)
  # The function should handle continuous variables gracefully, even though PMF is typically used for discrete variables.
  expect_equal(result, c(2/4, 1/4, 1/4), ignore_attr = TRUE)
  expect_equal(names(result), get_Y_t_names(Y_t_continuous))
})

test_that("get_pmf handles empty or NULL y_t", {
  Y_t <- list(values = c(0.1, 0.2, 0.3), type = 'continuous')
  result <- get_pmf(NULL, Y_t)
  expect_equal(result, rep(NaN, length(Y_t$values)), ignore_attr = TRUE)
  expect_equal(names(result), get_Y_t_names(Y_t))
})
