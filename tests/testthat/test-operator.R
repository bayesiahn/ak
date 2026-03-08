test_that("get_midpoints", {
  test_that("get_midpoints calculates midpoints correctly", {
    result <- get_midpoints(c(0.1, 0.3, 0.5))
    expect_equal(result, c(0.2, 0.4))
  })

  test_that("get_midpoints returns error for non-numeric input", {
    expect_error(get_midpoints(c("a", "b", "c")))
  })

  test_that("get_midpoints returns error for input of length less than 2", {
    expect_error(get_midpoints(c(0.1)))
  })
})

test_that("get_expectation_operator", {
  test_that("get_expectation_operator handles invalid Y_t", {
    expect_error(get_expectation_operator(NULL))
    expect_error(get_expectation_operator(list(values = c(1, 2, 3))))
    expect_error(get_expectation_operator(list(type = 'discrete')))
  })

  test_that("get_expectation_operator returns h_t for discrete type", {
    result <- get_expectation_operator(list(values = c(1, 2, 3), type = 'discrete'))
    expect_equal(result, c(1, 2, 3))
  })

  test_that("get_expectation_operator returns h_t for discrete type", {
    Y_t <- list(values = c(0.1, 0.3, 0.5), type = 'discrete',
                lower_bound = 0)
    result <- get_expectation_operator(Y_t)
    expect_equal(result, c(0.1, 0.3, 0.5))
  })

  test_that("get_expectation_operator returns h_t for continuous type", {
    Y_t <- list(values = c(0.1, 0.3, 0.5), type = 'continuous',
                lower_bound = 0)
    result <- get_expectation_operator(Y_t)
    expect_equal(result, c(0.05, 0.2, 0.4))
  })

  test_that("get_expectation_operator throws error for unknown type", {
    expect_error(get_expectation_operator(list(values = c(1, 2, 3), type = 'unknown')))
  })
})
