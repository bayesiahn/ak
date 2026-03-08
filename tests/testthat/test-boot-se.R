
test_that("models_to_standard_errors computes scalar standard errors", {
  set.seed(1234)
  y <- matrix(c(1, 2, 1, 1,
                1, 1, 1, 2,
                1, 1, 2, 2,
                1, 1, 1, 1), nrow = 4)
  g <- c(0, 3, 3, 3)
  opts <- list(bootstrap_counts = 20, mc_cores = 1)
  J <- 1

  models <- bootstrap_transition_models_from_wide_matrix(J, y, g, opts = opts)

  result <- models_to_standard_errors(models)

  # Expected standard errors (values updated after fix for zero-sum rows using uniform dist)
  expected_ATT_se <- c(0.1569428, 0.1101846)
  expect_equal(result$ATT, expected_ATT_se, tolerance = 1e-6)

  # Pre-treatment periods
  # Standard errors: 1 -> 1 in the first period has variation
  expect_true(result$Ps_control[[1]][[1]][1,1] != 0)
  expect_true(result$Ps_control[[1]][[1]][1,2] != 0)

  # Standard errors: 2 -> 2 in the first period has no variation (no obs.)
  expect_true(result$Ps_control[[1]][[1]][2,1] == 0)
  expect_true(result$Ps_control[[1]][[1]][2,2] == 0)

  # Post-treatment periods
  # Standard errors: 1 -> 1 in the second period has variation for treated
  expect_true(result$Ps_treated[[1]][[2]][1,1] != 0)
  expect_true(result$Ps_treated[[1]][[2]][1,2] != 0)

  # Standard errors: 2 -> 2 in the second period has no variation for treated
  expect_true(result$Ps_treated[[1]][[2]][2,1] == 0)
  expect_true(result$Ps_treated[[1]][[2]][2,2] == 0)

  # Standard errors: 1 -> 1 in the second period has no variation for control
  expect_true(result$Ps_control[[1]][[2]][1,1] == 0)
  expect_true(result$Ps_control[[1]][[2]][1,2] == 0)

  # Standard errors: 2 -> 2 in the second period has no variation for control (single obs.)
  expect_true(result$Ps_control[[1]][[2]][2,1] == 0)
  expect_true(result$Ps_control[[1]][[2]][2,2] == 0)
})


test_that("models_to_standard_errors computes scalar standard errors", {
  set.seed(1234)
  y <- matrix(c(1, 2, 1, 1,
                1, 1, 1, 2,
                1, 1, 2, 2,
                1, 1, 1, 1), nrow = 4)
  g <- c(0, 3, 3, 3)
  # Use explicit multistart_bootstrap_counts for reproducibility (test was written with this value)
  # Disable two-stage multistart to maintain backward compatibility with expected values
  opts <- list(bootstrap_counts = 20, mc_cores = 1, multistart_bootstrap_counts = 200,
               two_stage_multistart = FALSE)
  J <- 2

  models <- bootstrap_transition_models_from_wide_matrix(J, y, g, opts = opts)

  result <- models_to_standard_errors(models)

  # Expected standard errors (values updated after fix for zero-sum rows using uniform dist)
  expected_ATT <- c(0.2299218, 0.1259277)
  expect_equal(result$ATT, expected_ATT, tolerance = 1e-6)
})

