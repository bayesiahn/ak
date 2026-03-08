# =============================================================================
# Tests for Nonlinear Solver (first-stage-solver.R)
# =============================================================================

# --- Reparameterization Tests ---

test_that("simplex_to_unconstrained and unconstrained_to_simplex are inverses", {
  # Test with various probability vectors
  test_pmfs <- list(
    c(0.5, 0.5),
    c(0.3, 0.7),
    c(0.1, 0.2, 0.7),
    c(0.25, 0.25, 0.25, 0.25),
    c(0.01, 0.01, 0.98)
  )

  for (p in test_pmfs) {
    eta <- simplex_to_unconstrained(p)
    p_recovered <- unconstrained_to_simplex(eta)
    expect_equal(p_recovered, p, tolerance = 1e-10,
                 info = paste("Roundtrip failed for p =", paste(p, collapse = ", ")))
    expect_equal(length(eta), length(p) - 1)
  }
})

test_that("simplex_to_unconstrained handles edge cases", {
  # Single element (degenerate simplex)
  expect_equal(simplex_to_unconstrained(1), numeric(0))

  # Very small probabilities
  p_small <- c(1e-15, 1 - 1e-15)
  eta <- simplex_to_unconstrained(p_small)
  p_recovered <- unconstrained_to_simplex(eta)
  expect_equal(sum(p_recovered), 1)
  expect_true(all(p_recovered > 0))
})

test_that("unconstrained_to_simplex handles edge cases", {
  # Empty eta (single category)
  expect_equal(unconstrained_to_simplex(numeric(0)), 1)

  # Very large eta values (numerical stability)
  eta_large <- c(100, -100)
  p <- unconstrained_to_simplex(eta_large)
  expect_equal(sum(p), 1, tolerance = 1e-10)
  expect_true(all(p >= 0))
  expect_true(p[1] > 0.99)  # First element should dominate

  # All zeros -> uniform
  eta_zero <- c(0, 0)
  p <- unconstrained_to_simplex(eta_zero)
  expect_equal(p, rep(1/3, 3), tolerance = 1e-10)
})

# --- Vectorization/Unvectorization Roundtrip Tests ---

test_that("vectorize/unvectorize first-stage unconstrained is identity", {
  set.seed(1234)
  J <- 2
  K <- 2
  T_max <- 4
  treated_period <- 3

  # Generate a transition model
  transition_model <- generate_transition_model(
    outcomes = 1:K, J = J, T_max = T_max, treated_period = treated_period
  )
  sample <- generate_sample(transition_model, 100)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  # Run a few EM iterations to get a valid model
  result <- first_stage_from_guess(transition_model, data_for_est,
                                    opts = list(max_iter = 10))
  tm <- result$transition_model

  # Vectorize
  vec_result <- vectorize_first_stage_unconstrained(tm, data_for_est)
  theta <- vec_result$theta
  structure <- vec_result$structure

  # Unvectorize
  recovered <- unvectorize_first_stage_unconstrained(theta, structure)

  # Check priors roundtrip
  expect_equal(recovered$priors, tm$priors, tolerance = 1e-10)

  # Check p_y1d roundtrip
  for (j in 1:J) {
    expect_equal(recovered$p_y1d[[j]], tm$p_y1d[[j]], tolerance = 1e-10,
                 info = paste("p_y1d mismatch for group", j))
  }

  # Check Ps_control roundtrip
  for (j in 1:J) {
    for (t in seq_along(tm$Ps_control[[j]])) {
      expect_equal(recovered$Ps_control[[j]][[t]], tm$Ps_control[[j]][[t]],
                   tolerance = 1e-10,
                   info = paste("Ps_control mismatch for group", j, "period", t))
    }
  }

  # Check Ps_treated roundtrip
  for (j in 1:J) {
    for (t in seq_along(tm$Ps_treated[[j]])) {
      expect_equal(recovered$Ps_treated[[j]][[t]], tm$Ps_treated[[j]][[t]],
                   tolerance = 1e-10,
                   info = paste("Ps_treated mismatch for group", j, "period", t))
    }
  }
})

test_that("vectorize/unvectorize roundtrip works with K=3", {
  set.seed(42)
  J <- 2
  K <- 3
  T_max <- 4
  treated_period <- 3

  transition_model <- generate_transition_model(
    outcomes = 1:K, J = J, T_max = T_max, treated_period = treated_period
  )
  sample <- generate_sample(transition_model, 100)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  result <- first_stage_from_guess(transition_model, data_for_est,
                                    opts = list(max_iter = 10))
  tm <- result$transition_model

  vec_result <- vectorize_first_stage_unconstrained(tm, data_for_est)
  recovered <- unvectorize_first_stage_unconstrained(vec_result$theta, vec_result$structure)

  expect_equal(recovered$priors, tm$priors, tolerance = 1e-10)
  for (j in 1:J) {
    expect_equal(recovered$p_y1d[[j]], tm$p_y1d[[j]], tolerance = 1e-10)
    for (t in seq_along(tm$Ps_control[[j]])) {
      expect_equal(recovered$Ps_control[[j]][[t]], tm$Ps_control[[j]][[t]], tolerance = 1e-10)
      expect_equal(recovered$Ps_treated[[j]][[t]], tm$Ps_treated[[j]][[t]], tolerance = 1e-10)
    }
  }
})

# --- Objective Function Tests ---

test_that("first_stage_neg_log_likelihood matches EM log-likelihood", {
  set.seed(1234)
  J <- 2
  transition_model <- generate_transition_model(c(1, 2), J, 4, treated_period = 3)
  sample <- generate_sample(transition_model, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  # Run EM to convergence
  result <- first_stage_from_guess(transition_model, data_for_est,
                                    opts = list(max_iter = 50),
                                    compute_empirical_Ps = FALSE)
  tm <- result$transition_model
  em_ll <- tm$log_likelihood

  # Compute neg_ll at same point
  unit_weights <- rep(1, nrow(sample$y))
  vec_result <- vectorize_first_stage_unconstrained(tm, data_for_est)
  neg_ll <- first_stage_neg_log_likelihood(vec_result$theta, vec_result$structure,
                                            data_for_est, unit_weights)

  expect_equal(-neg_ll, em_ll, tolerance = 1e-6,
               info = "Solver objective should match EM log-likelihood at same parameter values")
})

# --- Solver Tests ---

test_that("solver achieves at least as good likelihood as EM starting point", {
  set.seed(1234)
  J <- 2
  transition_model <- generate_transition_model(c(1, 2), J, 4, treated_period = 3)
  sample <- generate_sample(transition_model, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  # Run a short EM (like the short-run stage)
  em_result <- first_stage_from_guess(transition_model, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)
  em_ll <- em_result$transition_model$log_likelihood

  # Run solver from EM output
  solver_result <- first_stage_from_solver(em_result$transition_model, data_for_est,
                                            opts = list(solver_max_iter = 200))
  solver_ll <- solver_result$transition_model$log_likelihood

  expect_true(solver_ll >= em_ll - 1e-6,
              info = paste("Solver LL =", solver_ll, "should be >= EM LL =", em_ll))
})

test_that("first_stage_from_solver returns correct output format", {
  set.seed(1234)
  J <- 2
  transition_model <- generate_transition_model(c(1, 2), J, 4, treated_period = 3)
  sample <- generate_sample(transition_model, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  em_result <- first_stage_from_guess(transition_model, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  solver_result <- first_stage_from_solver(em_result$transition_model, data_for_est,
                                            opts = list(solver_max_iter = 100))

  # Check output structure matches first_stage_from_guess output
  expect_true(is.list(solver_result))
  expect_true("transition_model" %in% names(solver_result))
  expect_true("posteriors" %in% names(solver_result))

  tm <- solver_result$transition_model
  expect_true("priors" %in% names(tm))
  expect_true("priors_treated" %in% names(tm))
  expect_true("p_y1d" %in% names(tm))
  expect_true("Ps_control" %in% names(tm))
  expect_true("Ps_treated" %in% names(tm))
  expect_true("log_likelihood" %in% names(tm))
  expect_true("Ps_control_from_pre" %in% names(tm))
  expect_true("Ps_treated_from_pre" %in% names(tm))
  expect_true("Ps_control_from_0" %in% names(tm))
  expect_true("Ps_treated_from_0" %in% names(tm))
  expect_true("pmfs_initial_control" %in% names(tm))
  expect_true("pmfs_initial_treated" %in% names(tm))

  # Check posteriors dimensions
  expect_equal(nrow(solver_result$posteriors), nrow(sample$y))
  expect_equal(ncol(solver_result$posteriors), J)

  # Check priors are valid probabilities
  expect_true(all(tm$priors >= 0))
  expect_equal(sum(tm$priors), 1, tolerance = 1e-10)
  expect_true(all(tm$priors_treated >= 0))
  expect_equal(sum(tm$priors_treated), 1, tolerance = 1e-10)

  # Check transition matrices are stochastic
  for (j in 1:J) {
    for (t in seq_along(tm$Ps_control[[j]])) {
      expect_equal(rowSums(tm$Ps_control[[j]][[t]]), rep(1, nrow(tm$Ps_control[[j]][[t]])),
                   tolerance = 1e-10)
    }
  }
})

# --- Integration Tests ---

test_that("first_stage with long_run_method='solver' works end-to-end", {
  set.seed(1234)
  J <- 2
  transition_model <- generate_transition_model(c(1, 2), J, 4, treated_period = 3)
  sample <- generate_sample(transition_model, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  opts_solver <- list(
    multistart_short_run_counts = 4,
    multistart_long_run_counts = 2,
    max_iter_short_run = 5,
    long_run_method = "solver",
    solver_max_iter = 100,
    two_stage_multistart = TRUE,
    parallel_multistart = FALSE
  )

  result <- first_stage(data_for_est, sample$g, J, opts = opts_solver)

  expect_true(!is.null(result$transition_model))
  expect_true(!is.null(result$posteriors))
  expect_true(!is.null(result$transition_model$log_likelihood))
  expect_true(result$transition_model$log_likelihood > -Inf)
})

test_that("solver produces comparable results to EM", {
  set.seed(1234)
  J <- 2
  transition_model <- generate_transition_model(c(1, 2), J, 4, treated_period = 3)
  sample <- generate_sample(transition_model, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  shared_opts <- list(
    multistart_short_run_counts = 6,
    multistart_long_run_counts = 2,
    max_iter_short_run = 10,
    two_stage_multistart = TRUE,
    parallel_multistart = FALSE
  )

  # EM long run
  opts_em <- c(shared_opts, list(long_run_method = "em", max_iter_long_run = 100))
  result_em <- first_stage(data_for_est, sample$g, J, opts = opts_em)

  # Solver long run
  opts_solver <- c(shared_opts, list(long_run_method = "solver", solver_max_iter = 200))
  result_solver <- first_stage(data_for_est, sample$g, J, opts = opts_solver)

  # Both should achieve reasonable log-likelihoods
  expect_true(result_em$transition_model$log_likelihood > -Inf)
  expect_true(result_solver$transition_model$log_likelihood > -Inf)

  # Solver should be at least as good (or very close) to EM
  # Allow small tolerance since multistart may pick different candidates
  ll_diff <- result_solver$transition_model$log_likelihood -
    result_em$transition_model$log_likelihood
  expect_true(ll_diff > -5,
              info = paste("Solver LL =", result_solver$transition_model$log_likelihood,
                           "EM LL =", result_em$transition_model$log_likelihood))
})

# --- Option Getter Tests ---

test_that("get_first_stage_long_run_method returns correct values", {
  expect_equal(get_first_stage_long_run_method(NULL), "solver")
  expect_equal(get_first_stage_long_run_method(list()), "solver")
  expect_equal(get_first_stage_long_run_method(list(long_run_method = "solver")), "solver")
  expect_equal(get_first_stage_long_run_method(list(long_run_method = "em")), "em")
})

test_that("get_solver_max_iter returns correct values", {
  expect_equal(get_solver_max_iter(NULL), 500)
  expect_equal(get_solver_max_iter(list(solver_max_iter = 200)), 200)
})

test_that("get_solver_rel_tol returns correct values", {
  expect_equal(get_solver_rel_tol(NULL), 1e-6)
  expect_equal(get_solver_rel_tol(list(solver_rel_tol = 1e-6)), 1e-6)
})

# --- Edge Case: J=1 ---

test_that("solver works with J=1", {
  set.seed(1234)
  J <- 1
  transition_model <- generate_transition_model(c(1, 2), J, 4, treated_period = 3)
  sample <- generate_sample(transition_model, 100)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  em_result <- first_stage_from_guess(transition_model, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  solver_result <- first_stage_from_solver(em_result$transition_model, data_for_est,
                                            opts = list(solver_max_iter = 100))

  expect_true(!is.null(solver_result$transition_model))
  expect_equal(length(solver_result$transition_model$priors), 1)
  expect_equal(solver_result$transition_model$priors, 1)
})

# =============================================================================
# Tests for Analytic Gradient
# =============================================================================

# --- softmax_gradient unit tests ---

test_that("softmax_gradient returns zero for uniform counts with uniform distribution", {
  probs <- c(0.25, 0.25, 0.25, 0.25)
  counts <- c(25, 25, 25, 25)
  grad <- softmax_gradient(counts, probs)
  expect_equal(grad, rep(0, 3), tolerance = 1e-10)
})

test_that("softmax_gradient returns correct values for known cases", {
  # Two categories: all weight in first
  probs <- c(0.5, 0.5)
  counts <- c(100, 0)
  grad <- softmax_gradient(counts, probs)
  # grad = counts[1] - probs[1] * sum(counts) = 100 - 0.5 * 100 = 50
  expect_equal(grad, 50)

  # Three categories: skewed
  probs <- c(0.2, 0.3, 0.5)
  counts <- c(20, 30, 50)  # Matches probs * 100
  grad <- softmax_gradient(counts, probs)
  expect_equal(grad, c(0, 0), tolerance = 1e-10)
})

test_that("softmax_gradient handles K=1 edge case", {
  expect_equal(softmax_gradient(10, 1), numeric(0))
})

# --- Option getter test ---

test_that("get_solver_use_gradient returns correct values", {
  expect_equal(get_solver_use_gradient(NULL), TRUE)
  expect_equal(get_solver_use_gradient(list()), TRUE)
  expect_equal(get_solver_use_gradient(list(solver_use_gradient = FALSE)), FALSE)
  expect_equal(get_solver_use_gradient(list(solver_use_gradient = TRUE)), TRUE)
})

# --- Finite-difference gradient verification ---

test_that("analytic gradient matches finite-difference gradient for J=1 K=2", {
  set.seed(1234)
  J <- 1; K <- 2; T_max <- 4; treated_period <- 3
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  # Analytic gradient
  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  # Finite-difference gradient (central differences)
  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient for J=2 K=2", {
  set.seed(42)
  J <- 2; K <- 2; T_max <- 4; treated_period <- 3
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient for J=2 K=3", {
  set.seed(99)
  J <- 2; K <- 3; T_max <- 4; treated_period <- 3
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 300)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient for J=3 K=2 T=5", {
  set.seed(555)
  J <- 3; K <- 2; T_max <- 5; treated_period <- 3
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 300)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient for J=2 K=2 T=6 treated_period=4 (more pretreatment periods)", {
  set.seed(600)
  J <- 2; K <- 2; T_max <- 6; treated_period <- 4
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient for J=1 K=3 T=4", {
  set.seed(601)
  J <- 1; K <- 3; T_max <- 4; treated_period <- 3
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient for J=3 K=3 T=5", {
  set.seed(602)
  J <- 3; K <- 3; T_max <- 5; treated_period <- 3
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 400)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient for J=2 K=2 T=8 treated_period=6 (many pretreatment periods)", {
  set.seed(603)
  J <- 2; K <- 2; T_max <- 8; treated_period <- 6
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 300)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient for J=2 K=3 T=6 treated_period=5 (late treatment)", {
  set.seed(604)
  J <- 2; K <- 3; T_max <- 6; treated_period <- 5
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 300)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient for J=2 K=4 T=4 (larger state space)", {
  set.seed(605)
  J <- 2; K <- 4; T_max <- 4; treated_period <- 3
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  sample <- generate_sample(tm, 400)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

test_that("analytic gradient matches finite-difference gradient with non-uniform weights", {
  set.seed(606)
  J <- 2; K <- 2; T_max <- 4; treated_period <- 3
  tm <- generate_transition_model(1:K, J, T_max, treated_period = treated_period)
  N <- 200
  sample <- generate_sample(tm, N)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- runif(N, 0.5, 2.0)

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  theta <- vec$theta
  structure <- vec$structure

  grad_analytic <- first_stage_neg_log_likelihood_gradient(theta, structure, data_for_est, unit_weights)

  h <- 1e-6
  grad_fd <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta; theta_plus[i] <- theta_plus[i] + h
    theta_minus <- theta; theta_minus[i] <- theta_minus[i] - h
    f_plus <- first_stage_neg_log_likelihood(theta_plus, structure, data_for_est, unit_weights)
    f_minus <- first_stage_neg_log_likelihood(theta_minus, structure, data_for_est, unit_weights)
    grad_fd[i] <- (f_plus - f_minus) / (2 * h)
  }

  expect_equal(grad_analytic, grad_fd, tolerance = 1e-4,
               info = paste("Max diff:", max(abs(grad_analytic - grad_fd))))
})

# --- Gradient near zero at converged EM solution ---

test_that("gradient is approximately zero at converged EM solution", {
  set.seed(1234)
  J <- 2
  tm <- generate_transition_model(1:2, J, 4, treated_period = 3)
  sample <- generate_sample(tm, 500)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)
  unit_weights <- rep(1, nrow(sample$y))

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 200, tol = 1e-8),
                                       compute_empirical_Ps = FALSE)

  vec <- vectorize_first_stage_unconstrained(em_result$transition_model, data_for_est)
  grad <- first_stage_neg_log_likelihood_gradient(vec$theta, vec$structure, data_for_est, unit_weights)

  expect_true(max(abs(grad)) < 0.5,
              info = paste("Max |gradient| =", max(abs(grad)), "at EM-converged point"))
})

# --- Solver with gradient achieves comparable likelihood ---

test_that("solver with analytic gradient achieves comparable likelihood to solver without", {
  set.seed(1234)
  J <- 2
  tm <- generate_transition_model(1:2, J, 4, treated_period = 3)
  sample <- generate_sample(tm, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  # Solver WITHOUT gradient
  result_no_grad <- first_stage_from_solver(
    em_result$transition_model, data_for_est,
    opts = list(solver_max_iter = 100, solver_use_gradient = FALSE)
  )

  # Solver WITH gradient
  result_grad <- first_stage_from_solver(
    em_result$transition_model, data_for_est,
    opts = list(solver_max_iter = 100, solver_use_gradient = TRUE)
  )

  # Both should achieve comparable likelihoods
  ll_diff <- abs(result_grad$transition_model$log_likelihood -
                   result_no_grad$transition_model$log_likelihood)
  expect_true(ll_diff < 1,
              info = paste("LL with gradient:", result_grad$transition_model$log_likelihood,
                           "LL without:", result_no_grad$transition_model$log_likelihood))
})

# --- J=1 edge case with gradient ---

test_that("solver with gradient works for J=1", {
  set.seed(1234)
  J <- 1
  tm <- generate_transition_model(1:2, J, 4, treated_period = 3)
  sample <- generate_sample(tm, 100)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  solver_result <- first_stage_from_solver(
    em_result$transition_model, data_for_est,
    opts = list(solver_max_iter = 100, solver_use_gradient = TRUE)
  )

  expect_true(!is.null(solver_result$transition_model))
  expect_equal(length(solver_result$transition_model$priors), 1)
  expect_equal(solver_result$transition_model$priors, 1)
})

# --- Integration test with gradient enabled ---

test_that("first_stage with solver and gradient works end-to-end", {
  set.seed(1234)
  J <- 2
  tm <- generate_transition_model(c(1, 2), J, 4, treated_period = 3)
  sample <- generate_sample(tm, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  opts_solver <- list(
    multistart_short_run_counts = 4,
    multistart_long_run_counts = 2,
    max_iter_short_run = 5,
    long_run_method = "solver",
    solver_max_iter = 100,
    solver_use_gradient = TRUE,
    two_stage_multistart = TRUE,
    parallel_multistart = FALSE
  )

  result <- first_stage(data_for_est, sample$g, J, opts = opts_solver)

  expect_true(!is.null(result$transition_model))
  expect_true(!is.null(result$posteriors))
  expect_true(result$transition_model$log_likelihood > -Inf)
})

# --- Speed comparison (informational) ---

test_that("solver with gradient is faster than without (informational)", {
  skip_on_cran()
  set.seed(1234)
  J <- 2; K <- 3; T_max <- 5
  tm <- generate_transition_model(1:K, J, T_max, treated_period = 3)
  sample <- generate_sample(tm, 200)
  data_for_est <- split_data_for_estimation(sample$y, sample$g)

  em_result <- first_stage_from_guess(tm, data_for_est,
                                       opts = list(max_iter = 10),
                                       compute_empirical_Ps = FALSE)

  time_no_grad <- system.time({
    first_stage_from_solver(em_result$transition_model, data_for_est,
                            opts = list(solver_max_iter = 50, solver_use_gradient = FALSE))
  })

  time_grad <- system.time({
    first_stage_from_solver(em_result$transition_model, data_for_est,
                            opts = list(solver_max_iter = 50, solver_use_gradient = TRUE))
  })

  message("Without gradient: ", round(time_no_grad["elapsed"], 2), "s")
  message("With gradient: ", round(time_grad["elapsed"], 2), "s")
  message("Speedup: ", round(time_no_grad["elapsed"] / max(time_grad["elapsed"], 0.01), 1), "x")
})
