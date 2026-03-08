test_that("estimate_post_em_params returns correct post-treatment parameters", {
  set.seed(1234)
  T_max <- 5
  J_max <- 1
  treated_period <- 3

  for (J in 1:J_max) {
    for (treated_period in 3:5) {
      # Setup mock data
      transition_model <- generate_transition_model(outcomes = c(1, 2), J = J, T_max = T_max)
      transition_model$treated_period <- treated_period
      sample <- generate_sample(transition_model, N = 100)
      y_indices_matrix <- sample$y
      g <- sample$g
      data_for_est <- split_data_for_estimation(y_indices_matrix, g)


      # Call estimate_post_em_params and verify the result
      posteriors <- transitions_to_posteriors(y_indices_matrix,
                                        transition_model$pmfs_initial_control,
                                        transition_model$Ps_control,
                                        transition_model$priors)
      post_em_params <- estimate_post_em_params(posteriors, data_for_est)

      # Check that post_em_params is a list and contains expected components
      expect_type(post_em_params, "list")
      expect_true("pmfs_pre_treatment_treated" %in% names(post_em_params))
      expect_true("prob_treated" %in% names(post_em_params))
      expect_true("priors_treated" %in% names(post_em_params))

      # Check the lengths of components of post_em_params
      expect_equal(length(post_em_params$pmfs_pre_treatment_treated), J)
    }
  }
})

test_that("compute_prob_treated handles valid inputs correctly", {
  # Test data
  g <- c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0)
  unit_weights <- runif(10)

  # Expected output
  expected_output <- sum(unit_weights[g>0]) / sum(unit_weights)

  # Compute function output
  prob_treated <- compute_prob_treated(g, unit_weights)

  # Check results
  expect_equal(prob_treated, expected_output)
})

test_that("compute_prob_treated handles NULL unit_weights", {
  # Test data
  g <- c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0)

  # Expected output
  expected_output <- mean(g > 0)

  # Compute function output
  prob_treated <- compute_prob_treated(g)

  # Check results
  expect_equal(prob_treated, expected_output)
})

test_that("compute_prob_treated handles edge cases", {
  # Case with no treated units
  g <- rep(0, 10)
  unit_weights <- runif(10)
  expect_equal(compute_prob_treated(g, unit_weights), 0)

  # Case with all treated units
  g <- rep(1, 10)
  expect_equal(compute_prob_treated(g, unit_weights), 1)

  # Case with inconsistent lengths of g and unit_weights
  g <- c(1, 0, 1, 0, 1)
  unit_weights <- runif(10) # Length mismatch
  expect_error(compute_prob_treated(g, unit_weights), "Length")

  # Case with empty input
  g <- numeric(0)
  unit_weights <- numeric(0)
  expect_error(compute_prob_treated(g, unit_weights), "empty")
})

test_that("compute_priors_treated returns correct priors for treated groups", {
  posteriors <- cbind(c(0, 1), c(1, 0))
  g <- c(0, 1)
  expected_priors <- c(1, 0)
  result <- compute_priors_treated(posteriors, g)
  expect_equal(result, expected_priors)

  posteriors <- cbind(c(1, 0), c(0, 1))
  g <- c(0, 1)
  expected_priors <- c(0, 1)
  result <- compute_priors_treated(posteriors, g)
  expect_equal(result, expected_priors)

  posteriors <- cbind(c(0, 1), c(1, 0))
  g <- c(1, 1)
  expected_priors <- c(0.5, 0.5)
  result <- compute_priors_treated(posteriors, g)
  expect_equal(result, expected_priors)

  posteriors <- cbind(c(1, 0), c(0, 1))
  unit_weights <- c(1, 2)
  g <- c(1, 1)
  expected_priors <- c(1/3, 2/3)
  result <- compute_priors_treated(posteriors, g, unit_weights)
  expect_equal(result, expected_priors)

  posteriors <- as.matrix(c(1, 1))
  g <- c(1, 1)
  expected_priors <- 1
  result <- compute_priors_treated(posteriors, g)
  expect_equal(result, expected_priors)
})

test_that("aggregate_ATT_by_latent_type correctly aggregates LTATTs", {
  ATT_list <- list(c(0.1, 0.2), c(0.3, 0.4))
  priors_treated <- c(0.8, 0.2)
  expected_agg_ATTs <- c(0.1 * 0.8 + 0.3 * 0.2, 0.2 * 0.8 + 0.4 * 0.2)

  result <- aggregate_ATT_by_latent_type(ATT_list, priors_treated)

  expect_equal(result, expected_agg_ATTs)
})

test_that("aggregate_ATT_by_latent_type handles mismatched ATT lengths", {
  ATT_list <- list(c(0.1, 0.2), c(0.3, 0.4, 0.5))
  priors_treated <- c(0.5, 0.5)

  expect_error(aggregate_ATT_by_latent_type(ATT_list, priors_treated))
})

test_that("estimate_post_em_params returns Ps_from_0 with first transition matching Ps_empirical", {
  set.seed(1234)
  T_max <- 5
  treated_period <- 3

  for (J in 1:2) {
    # Setup mock data
    transition_model <- generate_transition_model(outcomes = c(1, 2), J = J, T_max = T_max)
    transition_model$treated_period <- treated_period
    sample <- generate_sample(transition_model, N = 100)
    y_indices_matrix <- sample$y
    g <- sample$g
    data_for_est <- split_data_for_estimation(y_indices_matrix, g)

    # Compute posteriors and post-EM parameters
    posteriors <- transitions_to_posteriors(y_indices_matrix,
                                      transition_model$pmfs_initial_control,
                                      transition_model$Ps_control,
                                      transition_model$priors)
    post_em_params <- estimate_post_em_params(posteriors, data_for_est)

    # Check that Ps_from_0 components exist
    expect_true("Ps_control_from_0" %in% names(post_em_params))
    expect_true("Ps_treated_from_0" %in% names(post_em_params))

    # Check that Ps_from_0 have correct length (J latent groups)
    expect_equal(length(post_em_params$Ps_control_from_0), J)
    expect_equal(length(post_em_params$Ps_treated_from_0), J)

    # Check that Ps_from_0 have correct number of time periods (T-1)
    expect_equal(length(post_em_params$Ps_control_from_0[[1]]), T_max - 1)
    expect_equal(length(post_em_params$Ps_treated_from_0[[1]]), T_max - 1)

    # KEY TEST: First transition of Ps_from_0 should match first transition of Ps_empirical
    # This verifies that the non-Markovian and Markovian matrices agree at t=1
    for (j in seq_len(J)) {
      expect_equal(
        post_em_params$Ps_control_from_0[[j]][[1]],
        post_em_params$Ps_control_empirical[[j]][[1]],
        tolerance = 1e-10,
        info = paste("Ps_control first transition mismatch for j =", j)
      )
      expect_equal(
        post_em_params$Ps_treated_from_0[[j]][[1]],
        post_em_params$Ps_treated_empirical[[j]][[1]],
        tolerance = 1e-10,
        info = paste("Ps_treated first transition mismatch for j =", j)
      )
    }
  }
})

