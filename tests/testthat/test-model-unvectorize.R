library(testthat)

test_that("unvectorize_transition_model correctly reconstructs the transition model", {
  model <- list(
    priors = c(0.3, 0.7),
    priors_treated = c(0.4, 0.6),
    pmfs_initial_control = list(c(0.5, 0.5), c(0.2, 0.8)),
    pmfs_initial_treated = list(c(0.6, 0.4), c(0.3, 0.7)),
    pmfs_pre_treatment_control = list(c(0.7, 0.3), c(0.5, 0.5)),
    pmfs_pre_treatment_treated = list(c(0.8, 0.2), c(0.4, 0.6)),
    Ps_control_empirical = list(list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2, byrow = T)),
                                list(matrix(c(0.3, 0.7, 0.2, 0.8), nrow = 2, byrow = T))),
    Ps_treated_empirical = list(list(matrix(c(0.1, 0.9, 0.4, 0.6), nrow = 2, byrow = T)),
                                list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2, byrow = T)))
  )

  vectorized_model <- vectorize_transition_model(model)
  reconstructed_model <- unvectorize_transition_model(vectorized_model, model)

  expect_equal(reconstructed_model, model)
})

test_that("unvectorize_transition_model handles minimal = FALSE", {
  model_structure <- list(
    priors = c(0.3, 0.7),
    priors_treated = c(0.4, 0.6),
    pmfs_initial_control = list(c(0.5, 0.5)),
    pmfs_initial_treated = list(c(0.6, 0.4)),
    pmfs_pre_treatment_control = list(c(0.7, 0.3)),
    pmfs_pre_treatment_treated = list(c(0.8, 0.2)),
    Ps_control_empirical = list(list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2)),
                                list(matrix(c(0.6, 0.4, 0.3, 0.7), nrow = 2))),
    Ps_treated_empirical = list(list(matrix(c(0.6, 0.4, 0.3, 0.7), nrow = 2)),
                                list(matrix(c(0.3, 0.7, 0.3, 0.7), nrow = 2)))
  )

  vectorized_model <- vectorize_transition_model(model_structure, minimal = FALSE)
  reconstructed_model <- unvectorize_transition_model(vectorized_model, model_structure, minimal = FALSE)

  expect_equal(reconstructed_model, model_structure)
})

