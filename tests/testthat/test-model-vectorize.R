library(testthat)

test_that("vectorize_transition_model correctly vectorizes transition model components", {
    model <- list(
      priors = c(0.3, 0.7),
      priors_treated = c(0.4, 0.6),
      pmfs_initial_control = list(c(0.5, 0.5), c(0.2, 0.8)),
      pmfs_initial_treated = list(c(0.6, 0.4), c(0.3, 0.7)),
      pmfs_pre_treatment_control = list(c(0.7, 0.3), c(0.5, 0.5)),
      pmfs_pre_treatment_treated = list(c(0.8, 0.2), c(0.4, 0.6)),
      Ps_control_empirical = list(list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2))),
      Ps_treated_empirical = list(list(matrix(c(0.6, 0.4, 0.3, 0.7), nrow = 2)))
    )

    result <- vectorize_transition_model(model)

    expect_type(result, "list")
    expect_named(result, c("priors", "priors_treated", "pmfs_initial_control",
                           "pmfs_initial_treated", "pmfs_pre_treatment_control",
                           "pmfs_pre_treatment_treated", "Ps_control_empirical",
                           "Ps_treated_empirical"))

    expect_equal(result$priors, vectorize_pmf(model$priors))
    expect_equal(result$priors_treated, vectorize_pmf(model$priors_treated))
    expect_equal(result$pmfs_initial_control, vectorize_pmfs(model$pmfs_initial_control))
    expect_equal(result$pmfs_initial_treated, vectorize_pmfs(model$pmfs_initial_treated))
    expect_equal(result$pmfs_pre_treatment_control, vectorize_pmfs(model$pmfs_pre_treatment_control))
    expect_equal(result$pmfs_pre_treatment_treated, vectorize_pmfs(model$pmfs_pre_treatment_treated))
    expect_equal(result$Ps_control_empirical, vectorize_transition_matrices(model$Ps_control_empirical[[1]]))
    expect_equal(result$Ps_treated_empirical, vectorize_transition_matrices(model$Ps_treated_empirical[[1]]))
  })

  test_that("vectorize_transition_model handles minimal = FALSE correctly", {
    model <- list(
      priors = c(0.3, 0.7),
      priors_treated = c(0.4, 0.6),
      pmfs_initial_control = list(c(0.5, 0.5)),
      pmfs_initial_treated = list(c(0.6, 0.4)),
      pmfs_pre_treatment_control = list(c(0.7, 0.3)),
      pmfs_pre_treatment_treated = list(c(0.8, 0.2)),
      Ps_control_empirical = list(list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2)),
                                  list(matrix(c(0.6, 0.4, 0.3, 0.7), nrow = 2))),
      Ps_treated_empirical = list(list(matrix(c(0.6, 0.4, 0.3, 0.7), nrow = 2)),
                                  list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2)))
    )

    result <- vectorize_transition_model(model, minimal = FALSE)

    expect_equal(result$priors, vectorize_pmf(model$priors, minimal = FALSE))
    expect_equal(result$priors_treated, vectorize_pmf(model$priors_treated, minimal = FALSE))
    expect_equal((result$pmfs_initial_control), vectorize_pmfs(model$pmfs_initial_control, minimal = FALSE))
    expect_equal((result$pmfs_initial_treated), vectorize_pmfs(model$pmfs_initial_treated, minimal = FALSE))
    expect_equal((result$pmfs_pre_treatment_control), vectorize_pmfs(model$pmfs_pre_treatment_control, minimal = FALSE))
    expect_equal((result$pmfs_pre_treatment_treated), vectorize_pmfs(model$pmfs_pre_treatment_treated, minimal = FALSE))
    expect_equal((result$Ps_control_empirical),
                 vectorize_transition_matrices_list(model$Ps_control_empirical, minimal = FALSE))
    expect_equal((result$Ps_treated_empirical),
                 vectorize_transition_matrices_list(model$Ps_treated_empirical, minimal = FALSE))
  })

  test_that("vectorize_transition_model handles empty transition model", {
    model <- list(
      priors = numeric(0),
      priors_treated = numeric(0),
      pmfs_initial_control = list(),
      pmfs_initial_treated = list(),
      pmfs_pre_treatment_control = list(),
      pmfs_pre_treatment_treated = list(),
      Ps_control_empirical = list(),
      Ps_treated_empirical = list()
    )

    result <- vectorize_transition_model(model)

    expect_equal(result$priors, numeric(0))
    expect_equal(result$priors_treated, numeric(0))
    expect_equal(result$pmfs_initial_control, list())
    expect_equal(result$pmfs_initial_treated, list())
    expect_equal(result$pmfs_pre_treatment_control, list())
    expect_equal(result$pmfs_pre_treatment_treated, list())
    expect_equal(result$Ps_control_empirical, list())
    expect_equal(result$Ps_treated_empirical, list())
  })

  test_that("vectorize_transition_model handles single-column transition matrices", {
    model <- list(
      Ps_control_empirical = list(list(matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2))),
      Ps_treated_empirical = list(list(matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)))
    )

    result <- vectorize_transition_model(model)

    expect_equal(result$Ps_control_empirical, vectorize_transition_matrices(model$Ps_control_empirical[[1]]))
    expect_equal(result$Ps_treated_empirical, vectorize_transition_matrices(model$Ps_treated_empirical[[1]]))
})

test_that("completely_vectorize_transition_model correctly vectorizes a transition model", {
  model <- list(
    priors = c(0.3, 0.7),
    priors_treated = c(0.4, 0.6),
    pmfs_initial_control = list(c(0.5, 0.5), c(0.2, 0.8)),
    pmfs_initial_treated = list(c(0.6, 0.4), c(0.3, 0.7)),
    pmfs_pre_treatment_control = list(c(0.7, 0.3), c(0.5, 0.5)),
    pmfs_pre_treatment_treated = list(c(0.8, 0.2), c(0.4, 0.6)),
    Ps_control_empirical = list(list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2)),
                                list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2))),
    Ps_treated_empirical = list(list(matrix(c(0.6, 0.4, 0.3, 0.7), nrow = 2)),
                                list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2)))
  )

  expected_vectorized <- c(
    vectorize_pmf(model$priors),
    vectorize_pmf(model$priors_treated),
    vectorize_pmfs(model$pmfs_initial_control),
    vectorize_pmfs(model$pmfs_initial_treated),
    vectorize_pmfs(model$pmfs_pre_treatment_control),
    vectorize_pmfs(model$pmfs_pre_treatment_treated),
    vectorize_transition_matrices_list(model$Ps_control_empirical),
    vectorize_transition_matrices_list(model$Ps_treated_empirical)
  )

  result <- completely_vectorize_transition_model(model)
  expect_equal(unname(result), expected_vectorized)
})

test_that("completely_vectorize_transition_model handles minimal = FALSE", {
  model <- list(
    priors = c(0.3, 0.7),
    priors_treated = c(0.4, 0.6),
    pmfs_initial_control = list(c(0.5, 0.5), c(0.2, 0.8)),
    pmfs_initial_treated = list(c(0.6, 0.4), c(0.3, 0.7)),
    pmfs_pre_treatment_control = list(c(0.7, 0.3), c(0.5, 0.5)),
    pmfs_pre_treatment_treated = list(c(0.8, 0.2), c(0.4, 0.6)),
    Ps_control_empirical = list(list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2)),
                                list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2))),
    Ps_treated_empirical = list(list(matrix(c(0.6, 0.4, 0.3, 0.7), nrow = 2)),
                                list(matrix(c(0.5, 0.5, 0.2, 0.8), nrow = 2)))
  )

  expected_vectorized <- c(
    vectorize_pmf(model$priors, minimal = FALSE),
    vectorize_pmf(model$priors_treated, minimal = FALSE),
    vectorize_pmfs(model$pmfs_initial_control, minimal = FALSE),
    vectorize_pmfs(model$pmfs_initial_treated, minimal = FALSE),
    vectorize_pmfs(model$pmfs_pre_treatment_control, minimal = FALSE),
    vectorize_pmfs(model$pmfs_pre_treatment_treated, minimal = FALSE),
    vectorize_transition_matrices_list(model$Ps_control_empirical, minimal = FALSE),
    vectorize_transition_matrices_list(model$Ps_treated_empirical, minimal = FALSE)
  )

  result <- completely_vectorize_transition_model(model, minimal = FALSE)
  expect_equal(unname(result), expected_vectorized)
})

test_that("completely_vectorize_transition_model handles empty model", {
  model <- list(
    priors = numeric(0),
    priors_treated = numeric(0),
    pmfs_initial_control = list(),
    pmfs_initial_treated = list(),
    pmfs_pre_treatment_control = list(),
    pmfs_pre_treatment_treated = list(),
    Ps_control_empirical = list(),
    Ps_treated_empirical = list()
  )

  result <- completely_vectorize_transition_model(model)
  expect_equal(result, numeric(0))
})
