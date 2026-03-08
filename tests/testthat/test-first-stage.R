test_that("first_stage_from_guess converges to a solution", {
  # Setup: Create a dummy transition model guess and y_indices_matrix
  set.seed(1234)
  transition_model_guess <- generate_transition_model(c(1,2), 2, 3)
  y_indices_matrix <- generate_sample(transition_model_guess, 50)$y
  g <- c(rep(0, nrow(y_indices_matrix) - 1), 3)

  # Exercise: Run the function
  data_for_est <- split_data_for_estimation(y_indices_matrix, g)
  transition_model <- first_stage_from_guess(transition_model_guess, data_for_est, g)$transition_model

  # Verify: Check if a transition model is returned with expected components
  expect_true(is.list(transition_model))
  expect_true("priors" %in% names(transition_model))
  expect_true("pmfs_initial_control" %in% names(transition_model))
  expect_true("Ps_control" %in% names(transition_model))
  expect_true("log_likelihood" %in% names(transition_model))

  # Check if log likelihood is not -Inf (indicating some convergence)
  expect_true(transition_model$log_likelihood > -Inf)

  # Check if log likelihood improved
  pmfs_initial_control <- transition_model_guess$pmfs_initial_control
  Ps_control <- transition_model_guess$Ps_control
  priors <- transition_model_guess$priors
  transition_model_guess_likelihood <- transitions_to_log_likelihood(y_indices_matrix,
    pmfs_initial_control, Ps_control, priors)
  expect_true(transition_model$log_likelihood > transition_model_guess_likelihood)
})

test_that("generate_first_stage_guess returns a valid transition model", {
  y_indices_matrix <- matrix(sample(1:3, 9, replace = TRUE), nrow = 3)
  J <- 2
  result <- generate_first_stage_guess(y_indices_matrix, J)

  # Verify that the result is a list and contains expected components
  expect_true(is.list(result))
  expect_true("priors" %in% names(result))
  expect_true("pmfs_initial_control" %in% names(result))
  expect_true("Ps_control" %in% names(result))

  # Verify the dimensions of components match expected values
  expect_equal(length(result$priors), J)
  expect_true(all(sapply(result$pmfs_initial_control, length) == length(unique(c(y_indices_matrix)))))
  expect_equal(length(result$Ps_control), ncol(y_indices_matrix) - 1)
})

test_that("generate_first_stage_guess handles different numbers of outcomes", {
  # Test with different numbers of unique outcomes in y_indices_matrix
  for (num_outcomes in 1:2) {
    y_indices_matrix <- matrix(sample(1:num_outcomes, 10, replace = TRUE), nrow = 2)
    # Ensure at least some observations are contained in each outcome
    y_indices_matrix[1:num_outcomes, 1] <- 1:num_outcomes
    J <- 2
    result <- generate_first_stage_guess(y_indices_matrix, J)

    expect_equal(length(result$pmfs_initial_control), J)
    expect_true(all(sapply(result$pmfs_initial_control, length) == num_outcomes))
  }
})

test_that("first_stage returns a valid model", {
  set.seed(1234)
  J <- 2
  transition_model <- generate_transition_model(c(1,2), 2, 4, treated_period = 3)
  sample <- generate_sample(transition_model, 200)
  y_indices_matrix <- sample$y
  g <- sample$g
  opts <- list(multistart_counts = 2)
  data_for_est <- split_data_for_estimation(y_indices_matrix, g)

  result <- first_stage(data_for_est, g, J, opts = opts)

  # Check if a valid transition model is returned
  expect_true(is.list(result))
  expect_true("transition_model" %in% names(result))
  expect_true("posteriors" %in% names(result))
  expect_true("log_likelihood" %in% names(result$transition_model))
  expect_true("priors" %in% names(result$transition_model))
  expect_true("pmfs_initial_control" %in% names(result$transition_model))
  expect_true("Ps_control" %in% names(result$transition_model))

  # Make sure Ps_control_from_pre are valid
  treated_period <- max(sample$g)
  for (j in 1:J) {
    transition_model <- result$transition_model
    Ps_control_from_pre_j <- transition_model$Ps_control_from_pre[[j]]
    Ps_treated_from_pre_j <- transition_model$Ps_treated_from_pre[[j]]
    Ps_control_j <- transition_model$Ps_control[[j]]
    Ps_treated_j <- transition_model$Ps_treated[[j]]

    # Should be identical on treatment
    # (equivalent to Markovian in the first post-treatment period)
    expect_equal(Ps_control_from_pre_j[[treated_period - 1]], Ps_control_j[[treated_period - 1]], tolerance = 1e-2)
    expect_equal(Ps_treated_from_pre_j[[treated_period - 1]], Ps_treated_j[[treated_period - 1]], tolerance = 1e-2)

    # Might not be identical just after treatment
    expect_false(all(Ps_control_from_pre_j[[treated_period]] == Ps_control_j[[treated_period]]))
    expect_false(all(Ps_treated_from_pre_j[[treated_period]] == Ps_treated_j[[treated_period]]))
  }
})

test_that("first_stage raises error if y_indices_matrix has no columns", {
  y_indices_matrix <- matrix(numeric(0), nrow = 2, ncol = 0)
  g <- c(0, 1)

  data_for_est <- expect_error(split_data_for_estimation(y_indices_matrix, g),
                               "y_indices must have at least two columns.")
})

test_that("first_stage works if y_indices_matrix has at least two columns", {
  y_indices_matrix <- matrix(c(1, 2, 1, 2, 1, 2), nrow = 2)
  g <- c(0, 3)
  data_for_est <- split_data_for_estimation(y_indices_matrix, g)
  result <- first_stage(data_for_est, g, J = 2)

  # Check that result is a list with a transition_model
  expect_true(is.list(result))
  expect_true("transition_model" %in% names(result))
})

test_that("M_step returns correct components", {
  set.seed(1234)
  posteriors <- matrix(c(0.1, 0.1, 0.9, 0.9, 0.9, 0.1), ncol = 2)
  y_indices_matrix <- matrix(sample(1:2, 9, replace = TRUE), nrow = 3)
  g <- c(rep(0, nrow(y_indices_matrix) - 1), ncol(y_indices_matrix))
  data_for_est <- split_data_for_estimation(y_indices_matrix, g)
  result <- M_step(posteriors, data_for_est)

  expect_true(is.list(result))
  expect_true("priors" %in% names(result))
  expect_true("pmfs_initial_control" %in% names(result))
  expect_true("Ps_control" %in% names(result))
  expect_true("log_likelihood" %in% names(result))
})

test_that("M_step_priors returns normalized priors", {
  posteriors <- matrix(c(0.2, 0.9, 0.8, 0.1), ncol = 2)
  result <- M_step_priors(posteriors)

  expect_true(is.numeric(result))
  expect_equal(length(result), ncol(posteriors))
  expect_true(all(result >= 0 & result <= 1))
  expect_equal(sum(result), 1)
})

test_that("M_step_priors handles edge cases", {
  # Test for edge case: Empty matrix
  posteriors <- matrix(numeric(0), ncol = 0)
  result <- M_step_priors(posteriors)
  expect_equal(length(result), 0)

  # Test for edge case: Single column matrix
  posteriors <- matrix(runif(4), ncol = 1)
  result <- M_step_priors(posteriors)
  expect_equal(length(result), 1)
  expect_equal(result, 1)  # Single group should have a prior of 1
})

test_that("M_step_initial_outcome_dists returns correct distributions", {
  # Case 1: J = 1
  posteriors <- matrix(rep(1,3), ncol = 1)
  initial_outcome_indicators <- matrix(c(1, 0, 0, 0, 1, 1), ncol = 2)
  expected_result <- list(c(1/3, 2/3))
  result <- M_step_initial_outcome_dists(posteriors, initial_outcome_indicators)

  expect_equal(result, expected_result)

  # Case 2: J = 2
  posteriors <- matrix(c(1, 0.5, 0, 0, 0.5, 1), ncol = 2)
  initial_outcome_indicators <- matrix(c(1, 1, 0, 0, 0, 1), ncol = 2)
  expected_result <- list(c(1, 0.), c(1/3, 2/3))
  result <- M_step_initial_outcome_dists(posteriors, initial_outcome_indicators)

  expect_equal(result, expected_result)
})

test_that("M_step_initial_outcome_dists_j returns normalized distribution", {
  posteriors_j <- c(0.5, 0.5)
  initial_outcome_indicators <- matrix(c(1, 0, 0, 1), ncol = 2)
  result <- M_step_initial_outcome_dists_j(posteriors_j, initial_outcome_indicators)

  expect_equal(length(result), ncol(initial_outcome_indicators))
  expect_true(all(result >= 0 & result <= 1))
  expect_equal(sum(result), 1)
})


test_that("M_step_Ps returns correct structure", {
  posteriors <- matrix(runif(6), nrow = 2)
  transition_indicators <- list(list(matrix(1:4, nrow = 2), matrix(5:8, nrow = 2)))
  result <- M_step_Ps(posteriors, transition_indicators)

  expect_true(is.list(result))
  expect_equal(length(result), ncol(posteriors))
  expect_true(all(sapply(result, is.list)))
})

test_that("M_step_Ps_j returns correct structure", {
  posteriors_j <- runif(2)
  transition_indicators <- list(list(matrix(1:4, nrow = 2), matrix(5:8, nrow = 2)))
  result <- M_step_Ps_j(posteriors_j, transition_indicators)

  expect_true(is.list(result))
  expect_equal(length(result), length(transition_indicators))
  # Add more checks for the structure of each matrix in the list
})

test_that("M_step_P_jt returns normalized matrix", {
  posteriors_j <- runif(2)
  transition_indicators_t <- list(matrix(1:4, nrow = 2))
  result <- M_step_P_jt(posteriors_j, transition_indicators_t)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(transition_indicators_t[[1]]))
  expect_equal(ncol(result), ncol(transition_indicators_t[[1]]))

  # Check if rows sum to 1 for normalization
  expect_equal(apply(result, 1, sum), rep(1, nrow(result)))
})

test_that("M_step_P_jt returns expected value", {
  posteriors_j <- c(1, 0.3, 0.7)
  transition_indicators_t <- list(matrix(rep(0, 4), nrow = 2),
                                  matrix(rep(0, 4), nrow = 2),
                                  matrix(rep(0, 4), nrow = 2))
  transition_indicators_t[[1]][1,2] <- 1
  transition_indicators_t[[2]][2,1] <- 1
  transition_indicators_t[[3]][2,2] <- 1
  result <- M_step_P_jt(posteriors_j, transition_indicators_t)

  # Check if rows sum to 1 for normalization
  expect_equal(apply(result, 1, sum), rep(1, nrow(result)))

  # Check if matched with expected values
  expect_equal(result[1,1], 0)
  expect_equal(result[1,2], 1)
  expect_equal(result[2,1], 0.3)
  expect_equal(result[2,2], 0.7)
})

test_that("multistart generates different initial guesses", {
  # Test that generate_first_stage_guess produces different guesses
  # when called multiple times (as in first_stage multistart)
  set.seed(1234)
  y_indices_matrix <- matrix(sample(1:2, 20, replace = TRUE), nrow = 4)
  J <- 2
  n_guesses <- 5

  # Generate multiple guesses like first_stage does
  guesses <- lapply(seq_len(n_guesses), function(x) {
    generate_first_stage_guess(y_indices_matrix, J, treated_period = 3)
  })

  # Extract priors from each guess for comparison
  all_priors <- lapply(guesses, function(g) g$priors)

  # Check that not all guesses have identical priors
  # (at least some should be different due to random initialization)
  unique_priors <- unique(all_priors)
  expect_true(length(unique_priors) > 1,
              info = "Multistart should generate different initial guesses")

  # Also check transition matrices are different
  all_Ps <- lapply(guesses, function(g) g$Ps_control[[1]][[1]])
  unique_Ps <- unique(all_Ps)
  expect_true(length(unique_Ps) > 1,
              info = "Multistart should generate different transition matrices")
})

test_that("get_first_stage_multistart_counts scales with J", {
 # J=1 should return 1 (no latent groups to estimate)
 expect_equal(get_first_stage_multistart_counts(NULL, J = 1), 1)

 # J=2 should return 2000 (base_count * (J-1) = 2000 * 1)
 expect_equal(get_first_stage_multistart_counts(NULL, J = 2), 2000)

 # J=3 should return 4000 (base_count * (J-1) = 2000 * 2, capped at 4000)
 expect_equal(get_first_stage_multistart_counts(NULL, J = 3), 4000)

 # Explicit multistart_counts should override
 opts <- list(multistart_counts = 50)
 expect_equal(get_first_stage_multistart_counts(opts, J = 3), 50)

 # Custom base_count should work, but capped at max_count (4000)
 opts <- list(multistart_counts_base = 100)
 expect_equal(get_first_stage_multistart_counts(opts, J = 3), 200)  # min(100 * 2, 4000) = 200
})

# ============================================================================
# Two-Stage Multistart Tests
# ============================================================================

test_that("get_first_stage_short_run_counts returns correct values", {
  # Default should fall back to multistart_counts
  expect_equal(get_first_stage_short_run_counts(NULL, J = 2), 2000)
  expect_equal(get_first_stage_short_run_counts(NULL, J = 3), 4000)

  # Explicit short_run_counts should override
  opts <- list(multistart_short_run_counts = 100)
  expect_equal(get_first_stage_short_run_counts(opts, J = 2), 100)

  # Should use multistart_counts if short_run not specified
  opts <- list(multistart_counts = 50)
  expect_equal(get_first_stage_short_run_counts(opts, J = 2), 50)
})

test_that("get_first_stage_long_run_counts returns correct values", {
  # Default is 20
  expect_equal(get_first_stage_long_run_counts(NULL), 20)

  # Explicit long_run_counts should override
  opts <- list(multistart_long_run_counts = 10)
  expect_equal(get_first_stage_long_run_counts(opts), 10)
})

test_that("get_first_stage_max_iter_short_run returns correct values", {
  # Default is 200
  expect_equal(get_first_stage_max_iter_short_run(NULL), 200)

  # Explicit value should override
  opts <- list(max_iter_short_run = 30)
  expect_equal(get_first_stage_max_iter_short_run(opts), 30)
})

test_that("get_first_stage_max_iter_long_run returns correct values", {
  # Default is 200
  expect_equal(get_first_stage_max_iter_long_run(NULL), 200)

  # Explicit value should override
  opts <- list(max_iter_long_run = 1000)
  expect_equal(get_first_stage_max_iter_long_run(opts), 1000)
})

test_that("use_two_stage_multistart returns correct values", {
  # Default is TRUE
  expect_true(use_two_stage_multistart(NULL))
  expect_true(use_two_stage_multistart(list()))

  # Explicit FALSE should disable
  opts <- list(two_stage_multistart = FALSE)
  expect_false(use_two_stage_multistart(opts))

  # Explicit TRUE should enable
  opts <- list(two_stage_multistart = TRUE)
  expect_true(use_two_stage_multistart(opts))
})

test_that("select_top_candidates selects correct number of models", {
  # Create mock estimated models
  mock_models <- lapply(1:10, function(i) {
    list(transition_model = list(log_likelihood = -i))  # Higher i = lower likelihood
  })

  # Select top 3
  top3 <- select_top_candidates(mock_models, 3)
  expect_equal(length(top3), 3)

  # Should be the ones with highest likelihood (lowest i, i.e., -1, -2, -3)
  likelihoods <- sapply(top3, function(x) x$transition_model$log_likelihood)
  expect_equal(sort(likelihoods, decreasing = TRUE), c(-1, -2, -3))

  # When k > length, should return all
  all_models <- select_top_candidates(mock_models, 100)
  expect_equal(length(all_models), 10)
})

test_that("two-stage and single-stage first_stage both work", {
  # Generate test data (same pattern as test "first_stage returns a valid model")
  set.seed(1234)
  J <- 2
  transition_model <- generate_transition_model(c(1, 2), J, 4, treated_period = 3)
  sample <- generate_sample(transition_model, 200)
  y_indices_matrix <- sample$y
  g <- sample$g

  data_for_est <- split_data_for_estimation(y_indices_matrix, g)

  # Test single-stage (original behavior)
  opts_single <- list(
    multistart_counts = 3,
    max_iter = 20,
    two_stage_multistart = FALSE,
    parallel_multistart = FALSE
  )
  result_single <- first_stage(data_for_est, g, J, opts = opts_single)

  expect_true(!is.null(result_single$transition_model))
  expect_true(!is.null(result_single$posteriors))
  expect_true(!is.null(result_single$transition_model$log_likelihood))

  # Test two-stage
  opts_two_stage <- list(
    multistart_short_run_counts = 6,
    multistart_long_run_counts = 2,
    max_iter_short_run = 5,
    max_iter_long_run = 20,
    two_stage_multistart = TRUE,
    parallel_multistart = FALSE
  )
  result_two_stage <- first_stage(data_for_est, g, J, opts = opts_two_stage)

  expect_true(!is.null(result_two_stage$transition_model))
  expect_true(!is.null(result_two_stage$posteriors))
  expect_true(!is.null(result_two_stage$transition_model$log_likelihood))
})

test_that("J=1 uses single-stage even when two_stage_multistart=TRUE", {
  # For J=1, there's no latent structure to estimate, so two-stage is not beneficial
  set.seed(1234)
  J <- 1
  transition_model <- generate_transition_model(c(1, 2), J, 4, treated_period = 3)
  sample <- generate_sample(transition_model, 100)
  y_indices_matrix <- sample$y
  g <- sample$g

  data_for_est <- split_data_for_estimation(y_indices_matrix, g)

  # Even with two_stage_multistart=TRUE, J=1 should use single-stage
  opts <- list(two_stage_multistart = TRUE, parallel_multistart = FALSE, max_iter = 20)
  result <- first_stage(data_for_est, g, J, opts = opts)

  expect_true(!is.null(result$transition_model))
  expect_equal(length(result$transition_model$priors), 1)
})
