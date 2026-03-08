test_that("difference_by_transitions returns correct difference", {
  # K = 1
  pmf <- 1
  P_treated <- matrix(c(1), nrow = 1)
  P_control <- matrix(c(1), nrow = 1)
  Y_values <- 1

  difference <- difference_by_transitions(pmf, P_treated, P_control, Y_values)

  # Calculate expected difference manually
  expected_difference <- as.matrix(0)
  expect_equal(difference, expected_difference)

  # K = 2
  pmf <- c(0.5, 0.5)
  P_treated <- matrix(c(0.3, 0.4, 0.7, 0.6), nrow = 2)
  P_control <- matrix(c(0.6, 0.5, 0.4, 0.5), nrow = 2)
  Y_values <- c(0, 1)

  difference <- difference_by_transitions(pmf, P_treated, P_control, Y_values)

  # Calculate expected difference manually
  expected_difference <- as.matrix(0.2)
  expect_equal(difference, expected_difference)
})

test_that("compute_LTATTs_under_markovian returns correct LTATTs", {
  pmf <- c(0.5, 0.5)
  Ps_treated <- list(matrix(c(0.7, 0.4, 0.3, 0.6), nrow = 2), matrix(c(0.6, 0.3, 0.4, 0.7), nrow = 2))
  Ps_control <- list(matrix(c(0.3, 0.4, 0.7, 0.6), nrow = 2), matrix(c(0.6, 0.5, 0.4, 0.5), nrow = 2))
  Y_values <- c(0, 1)

  LTATTs <- compute_LTATTs_under_markovian(pmf, Ps_treated, Ps_control, Y_values)

  # Verify the LTATTs are calculated correctly
  # Add specific checks here based on the expected behavior of your function
  expect_equal(LTATTs[1], -0.2)
  expect_equal(LTATTs[2], 0.07)
  expect_equal(length(LTATTs), length(Ps_treated))
})

test_that("compute_LTATTs_under_markovian handles mismatched lengths of Ps_treated and Ps_control", {
  pmf <- c(0.5, 0.5)
  Ps_treated <- list(matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2))  # Only one matrix
  Ps_control <- list(matrix(c(0.6, 0.4, 0.5, 0.5), nrow = 2), matrix(c(0.5, 0.5, 0.4, 0.6), nrow = 2))

  expect_error(compute_LTATTs_under_markovian(pmf, Ps_treated, Ps_control, Y_values))
})

test_that("get_expected_future_outcome computes correctly", {
  pmf <- c(0.5, 0.5)
  P_future <- matrix(c(0.1, -0.1, -0.1, 0.1), nrow = 2)
  Y_values <- c(0, 1)

  expected_outcome <- get_expected_future_outcome(pmf, P_future, Y_values)

  # Expected calculation: (0.5 * 0.1 + 0.5 * -0.1) * 0 + (0.5 * -0.1 + 0.5 * 0.1) * 1 = 0
  expect_equal(as.numeric(expected_outcome), 0)
})

test_that("get_expected_future_outcome handles different probabilities", {
  pmf <- c(0.3, 0.7)
  P_future <- matrix(c(0.2, 0.1, 0.4, 0.3), nrow = 2)
  Y_values <- c(1, 2)

  expected_outcome <- get_expected_future_outcome(pmf, P_future, Y_values)

  # Manual calculation for expected outcome
  manual_expected <- sum(pmf * (P_future %*% Y_values))

  expect_equal(as.numeric(expected_outcome), manual_expected)
})

test_that("get_expected_outcome computes correctly", {
  pmf <- c(0.5, 0.5)
  Y_values <- c(0, 1)

  expected_outcome <- get_expected_outcome(pmf, Y_values)

  # Expected calculation: (0.5 * 0) + (0.5 * 1) = 0.5
  expect_equal(as.numeric(expected_outcome), 0.5)
})

test_that("get_expected_outcome handles different pmfs correctly", {
  pmf <- c(0.3, 0.7)
  Y_values <- c(1, 3)

  expected_outcome <- get_expected_outcome(pmf, Y_values)

  # Manual calculation for expected outcome: (0.3 * 1) + (0.7 * 3) = 2.4
  expect_equal(as.numeric(expected_outcome), 2.4)
})


