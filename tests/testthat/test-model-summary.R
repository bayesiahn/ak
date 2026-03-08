library(testthat)

test_that("summarize_transition_model computes LTATTs and aggregates ATT correctly", {
  # Mock TransitionModel
  model <- generate_transition_model(J = 2, T_max = 3, treated_period = 3)
  result <- summarize_transition_model(model)

  # Expected LTATTs and LTATTs_from_0
  expected_LTATTs <- list(
    compute_LTATTs_under_markovian(model$pmfs_pre_treatment_treated[[1]],
                model$Ps_treated[[1]][2:2],
                model$Ps_control[[1]][2:2],
                model$Y$values),
    compute_LTATTs_under_markovian(model$pmfs_pre_treatment_treated[[2]],
                model$Ps_treated[[2]][2:2],
                model$Ps_control[[2]][2:2],
                model$Y$values)
  )
  expected_LTATTs_from_0 <- list(
    compute_LTATTs_under_markovian(model$pmfs_initial_treated[[1]],
                model$Ps_treated[[1]],
                model$Ps_control[[1]],
                model$Y$values),
    compute_LTATTs_under_markovian(model$pmfs_initial_treated[[2]],
                model$Ps_treated[[2]],
                model$Ps_control[[2]],
                model$Y$values)
  )

  # Expected aggregated ATT
  expected_ATT <- aggregate_ATT_by_latent_type(expected_LTATTs, model$priors_treated)
  expected_ATT_from_0 <- aggregate_ATT_by_latent_type(expected_LTATTs_from_0, model$priors_treated)

  # Verify LTATTs
  expect_equal(result$LTATTs, expected_LTATTs)

  # Verify aggregated ATT
  expect_equal(result$ATT, expected_ATT)

  # Verify that the model is returned
  expect_equal(result$model, model)
})

test_that("summarize_transition_model handles single latent group correctly", {
  # Mock TransitionModel with one group
  model <- generate_transition_model(J = 1, T_max = 2, treated_period = 2)

  result <- summarize_transition_model(model)

  # Expected LTATTs
  expected_LTATTs <- list(model$pmfs_pre_treatment_treated[[1]] %*%
    (model$Ps_treated[[1]][[1]] - model$Ps_control[[1]][[1]]) %*%
    c(0, 1) %>%
    as.numeric()
  )

  # Verify Ps_diffs
  Ps_diffs_expected <- model$Ps_treated[[1]][[1]] - model$Ps_control[[1]][[1]]
  expect_equal(result$Ps_diffs[[1]][[1]], Ps_diffs_expected)
  expect_equal(result$Ps_empirical_diffs[[1]][[1]], Ps_diffs_expected)

  # Expected aggregated ATT
  expected_ATT <- aggregate_ATT_by_latent_type(expected_LTATTs, model$priors_treated)

  # Verify LTATTs
  expect_equal(result$LTATTs, expected_LTATTs)

  # Verify aggregated ATT
  expect_equal(result$ATT, expected_ATT)

  # Verify that the model is returned
  expect_equal(result$model, model)
})
