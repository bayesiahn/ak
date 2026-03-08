set.seed(1234)

test_that("transition_to_log_likelihood_j returns correct results", {
  # Setup
  y_indices <- c(1, 2, 1)
  initial_outcome_dist_j <- c(0.5, 0.5)
  Ps_control_j <- list(matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2),
                                matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2))
  expected_result <- -2.079442

  # Exercise
  result <- transition_to_log_likelihood_j(y_indices,
    initial_outcome_dist_j, Ps_control_j)

  # Verify
  expect_equal(result, expected_result, tolerance = 1e-5)
})


test_that("transitions_to_log_likelihoods_j returns correct results", {
  set.seed(1234)

  # Setup
  y_indices_matrix <- matrix(c(1, 2, 1, 2, 2, 2), nrow = 2)
  initial_outcome_dist_j <- c(0.5, 0.5)
  Ps_control_j <- list(matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2),
   matrix(c(0.7, 0.5, 0.3, 0.5), nrow = 2))
  expected_result <- c(-2.590267, -2.079442)

  # Exercise
  result <- transitions_to_log_likelihoods_j(y_indices_matrix,
                                   initial_outcome_dist_j, Ps_control_j)
  result_by_each <- apply(y_indices_matrix, 1, function(y_indices) {
    transition_to_log_likelihood_j(y_indices, initial_outcome_dist_j, Ps_control_j)
  })

  # Verify
  expect_equal(result, expected_result, tolerance = 1e-5)
  expect_equal(result, result_by_each, tolerance = 1e-5)
})



test_that("transition_to_log_likelihood returns correct results", {
  set.seed(1234)

  # Setup
  y_indices <- c(1, 2, 1)

  # Setup: 1 latent group
  priors <- c(1)
  transition_model <- generate_transition_model(c(0,1), length(priors), length(y_indices),
                                    treated_period = 0)
  pmfs_initial_control <- transition_model$pmfs_initial_control
  Ps_control <- transition_model$Ps_control
  expected_result <- -2.120264

  # Exercise
  result <- transition_to_log_likelihood(y_indices,
    pmfs_initial_control, Ps_control, priors)
  single_group_result <- as.numeric(transition_to_log_likelihood_j(y_indices,
    pmfs_initial_control[[1]], Ps_control[[1]]))

  # Verify
  expect_equal(result, expected_result, tolerance = 1e-5)
  expect_equal(result, single_group_result, tolerance = 1e-5)
})

test_that("transitions_to_log_likelihood throws error for mismatched unit_weights", {
  y_indices_matrix <- matrix(c(1, 2, 1, 2), nrow = 2)
  priors <- c(0.3, 0.7)
  initial_outcome_dists <- list(rep(0.5, 2), rep(0.5, 2))
  Ps_control <- list(matrix(0.5, 2, 2), matrix(0.5, 2, 2))
  unit_weights <- c(1, 2, 3)  # Incorrect length

  expect_error(
    transitions_to_log_likelihood(y_indices_matrix, initial_outcome_dists, Ps_control, priors, unit_weights),
    "Length of unit_weights must match the number of rows in y_indices_matrix."
  )
})

test_that("transitions_to_log_likelihood works with no units", {
  y_indices_matrix <- matrix(numeric(0), nrow = 0, ncol = 3)
  priors <- c(0.3, 0.7)
  initial_outcome_dists <- list(rep(0.5, 2), rep(0.5, 2))
  Ps_control <- list(matrix(0.5, 2, 2), matrix(0.5, 2, 2))

  result <- transitions_to_log_likelihood(y_indices_matrix, initial_outcome_dists, Ps_control, priors)

  # Expect 0 log likelihood for empty data
  expect_equal(result, 0)
})

test_that("transitions_to_posteriors returns correct results", {
  # Setup
  set.seed(123456)
  y_indices_matrix <- rbind(c(1, 2, 1), c(2, 1, 2))
  priors <- c(0.7, 0.3)
  transition_model <- generate_transition_model(c(0,1), length(priors),
                                    ncol(y_indices_matrix),
                                    treated_period = 0)
  pmfs_initial_control <- transition_model$pmfs_initial_control
  Ps_control <- transition_model$Ps_control
  expected_result <- -4.04831

  # Exercise
  ll <- transitions_to_log_likelihood(y_indices_matrix,
                                pmfs_initial_control, Ps_control, priors)
  ll_each <- apply(y_indices_matrix, 1, function(y_indices) {
    transition_to_log_likelihood(y_indices,
                           pmfs_initial_control, Ps_control, priors)
  })

  # Verify
  expect_equal(ll, expected_result, tolerance = 1e-5)
  expect_equal(ll, sum(ll_each), tolerance = 1e-5)
})

test_that("transitions_to_posteriors returns correct results", {
  # Setup
  set.seed(1234)
  T_max <- 10
  N <- 200
  get_posteriors_from_generated_sample <- function(transition_model) {
    priors <- transition_model$priors
    pmfs_initial_control <- transition_model$pmfs_initial_control
    Ps_control <- transition_model$Ps_control
    sample <- generate_sample(transition_model, N)
    y_indices_matrix <- sample$y %>%
      y_to_y_indices(transition_model$Y$values)
    posteriors <- transitions_to_posteriors(y_indices_matrix, pmfs_initial_control,
                        Ps_control, priors)

    list(sample = sample, posteriors = posteriors)
  }

  # Setup: 1 latent group
  J <- 1
  transition_model <- generate_transition_model(c(0,1), J, T_max,
                                    treated_period = 0)
  out <- get_posteriors_from_generated_sample(transition_model)
  posteriors <- out$posteriors
  N_js <- out$sample$N_js

  # Verify if posteriors are proper posteriors
  ## Check dimension
  expect_equal(dim(posteriors), c(N, J))
  ## Check if each row is a vector of probabilities
  expect_equal(as.numeric(apply(posteriors, 1, sum)), rep(1, N))

  # Setup: 2 latent group
  J <- 2
  transition_model <- generate_transition_model(c(0,1), J, T_max,
                                    treated_period = 0)
  out <- get_posteriors_from_generated_sample(transition_model)
  posteriors <- out$posteriors
  N_js <- out$sample$N_js

  # Verify if posteriors are proper posteriors
  ## Check dimension
  expect_equal(dim(posteriors), c(N, J))
  ## Check if each row is a vector of probabilities
  expect_equal(as.numeric(apply(posteriors, 1, sum)), rep(1, N))

  # Check if posteriors provide meaningful info
  indices_1 <- 1:N_js[1]
  indices_2 <- (N_js[1]+1):sum(N_js)
  expect_true(mean((posteriors[,1] > posteriors[,2])[1:N_js[1]]) > 0.6)
  expect_true(mean((posteriors[,1] < posteriors[,2])[indices_2]) > 0.6)
})

