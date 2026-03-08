# =============================================================================
# Nonlinear Solver for First-Stage Estimation
# =============================================================================
#
# Provides an alternative to EM for the long-run stage of two-stage multistart.
# Uses stats::nlminb with softmax reparameterization to directly maximize
# the observed-data log-likelihood.

# --- Reparameterization: Probability Simplex <-> Unconstrained Reals ---

#' Convert Probability Simplex to Unconstrained Parameters
#'
#' Maps a K-dimensional probability vector (summing to 1) to K-1 unconstrained
#' real values using the multinomial logit (softmax inverse) transform:
#' \eqn{\eta_k = \log(p_k / p_K)} for \eqn{k = 1, \ldots, K-1}.
#'
#' @param p Numeric vector summing to 1 (a probability mass function).
#' @return Numeric vector of length \code{length(p) - 1}.
#' @export
simplex_to_unconstrained <- function(p) {
  K <- length(p)
  if (K <= 1) return(numeric(0))
  # Clamp to avoid log(0)
  p_clamped <- pmax(p, .Machine$double.eps)
  log(p_clamped[1:(K - 1)] / p_clamped[K])
}

#' Convert Unconstrained Parameters Back to Probability Simplex
#'
#' Inverse of \code{simplex_to_unconstrained}. Maps K-1 unconstrained real
#' values back to a K-dimensional probability vector using the softmax transform.
#' Uses the log-sum-exp trick for numerical stability.
#'
#' @param eta Numeric vector of length K-1.
#' @return Numeric vector of length K summing to 1.
#' @export
unconstrained_to_simplex <- function(eta) {
  if (length(eta) == 0) return(1)
  # Append 0 for the reference category (K-th element)
  eta_full <- c(eta, 0)
  # Log-sum-exp trick for numerical stability
  max_eta <- max(eta_full)
  exp_eta <- exp(eta_full - max_eta)
  exp_eta / sum(exp_eta)
}

# --- Softmax Gradient Helpers ---

#' Compute Gradient in Unconstrained (Softmax) Space
#'
#' Given weighted counts (sufficient statistics) and current simplex
#' probabilities, computes the gradient of the log-likelihood w.r.t.
#' the unconstrained softmax parameters using the identity:
#' \eqn{dL/d\eta_m = n_m - p_m \sum_k n_k} for \eqn{m = 1, \ldots, K-1}.
#'
#' @param weighted_counts Numeric vector of length K (weighted counts per category).
#' @param probs Numeric vector of length K (current simplex probabilities, sums to 1).
#' @return Numeric vector of length K-1 (gradient w.r.t. unconstrained params).
#' @keywords internal
#' @export
softmax_gradient <- function(weighted_counts, probs) {
  K <- length(probs)
  if (K <= 1) return(numeric(0))
  total <- sum(weighted_counts)
  (weighted_counts - probs * total)[1:(K - 1)]
}

#' Compute Weighted Transition Counts
#'
#' For a given time period, computes the K x K matrix of weighted
#' transition counts \eqn{n_{k,m} = \sum_i w_i \cdot 1\{Y_{i,t}=k, Y_{i,t+1}=m\}}.
#'
#' @param y_from Integer vector of origin states (length N).
#' @param y_to Integer vector of destination states (length N).
#' @param weights Numeric vector of weights (length N), typically \eqn{w_i \cdot Z_{ij}}.
#' @param K Integer number of outcome states.
#' @return K x K matrix of weighted transition counts.
#' @keywords internal
#' @export
compute_weighted_transition_counts <- function(y_from, y_to, weights, K) {
  counts <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    mask_k <- (y_from == k)
    if (any(mask_k)) {
      for (m in 1:K) {
        counts[k, m] <- sum(weights[mask_k & (y_to == m)])
      }
    }
  }
  counts
}

# --- First-Stage Parameter Vectorization (Unconstrained) ---

#' Vectorize First-Stage Parameters to Unconstrained Vector
#'
#' Converts all first-stage transition model parameters (priors, p_y1d,
#' transition matrices) into a single unconstrained real-valued vector
#' using softmax reparameterization.
#'
#' The parameter vector is assembled in this order:
#' \enumerate{
#'   \item \code{priors}: J-1 unconstrained values
#'   \item \code{p_y1d}: J blocks of (2K-1) values (joint P(Y_1, D | Z_j))
#'   \item \code{Ps_pre}: shared pretreatment transitions
#'   \item \code{Ps_control_post}: post-treatment control transitions
#'   \item \code{Ps_treated_post}: post-treatment treated transitions
#' }
#'
#' @param transition_model The transition model (from EM M-step output).
#' @param data_for_est Data structure from \code{split_data_for_estimation()}.
#' @return A list with:
#'   \describe{
#'     \item{theta}{Numeric vector of unconstrained parameters.}
#'     \item{structure}{Metadata needed for \code{unvectorize_first_stage_unconstrained()}.}
#'   }
#' @keywords internal
#' @export
vectorize_first_stage_unconstrained <- function(transition_model, data_for_est) {
  priors <- transition_model$priors
  p_y1d <- transition_model$p_y1d
  Ps_control <- transition_model$Ps_control
  Ps_treated <- transition_model$Ps_treated

  J <- length(priors)
  K <- nrow(p_y1d[[1]])  # Number of outcome states
  treated_period <- data_for_est$treated_period
  T_max <- ncol(data_for_est$y_indices_matrix)
  n_pre_periods <- max(treated_period - 2, 0)
  n_post_periods <- T_max - max(treated_period - 1, 1)

  # 1. Priors
  theta_priors <- simplex_to_unconstrained(priors)

  # 2. p_y1d: each p_y1d[[j]] is a K x 2 matrix; flatten column-major as a 2K-simplex
  theta_p_y1d <- unlist(lapply(p_y1d, function(p_y1d_j) {
    simplex_to_unconstrained(c(p_y1d_j))
  }))

  # 3. Ps_pre (shared pretreatment transitions, periods 1 to treated_period-2)
  theta_Ps_pre <- numeric(0)
  if (n_pre_periods > 0) {
    for (j in 1:J) {
      for (t in 1:n_pre_periods) {
        P_jt <- Ps_control[[j]][[t]]  # Same as Ps_treated[[j]][[t]] in pretreatment
        for (k in 1:nrow(P_jt)) {
          theta_Ps_pre <- c(theta_Ps_pre, simplex_to_unconstrained(P_jt[k, ]))
        }
      }
    }
  }

  # 4. Ps_control_post (post-treatment control transitions)
  theta_Ps_control_post <- numeric(0)
  if (n_post_periods > 0) {
    for (j in 1:J) {
      for (t_idx in 1:n_post_periods) {
        t <- n_pre_periods + t_idx  # Index into full Ps list
        P_jt <- Ps_control[[j]][[t]]
        for (k in 1:nrow(P_jt)) {
          theta_Ps_control_post <- c(theta_Ps_control_post, simplex_to_unconstrained(P_jt[k, ]))
        }
      }
    }
  }

  # 5. Ps_treated_post (post-treatment treated transitions)
  theta_Ps_treated_post <- numeric(0)
  if (n_post_periods > 0) {
    for (j in 1:J) {
      for (t_idx in 1:n_post_periods) {
        t <- n_pre_periods + t_idx  # Index into full Ps list
        P_jt <- Ps_treated[[j]][[t]]
        for (k in 1:nrow(P_jt)) {
          theta_Ps_treated_post <- c(theta_Ps_treated_post, simplex_to_unconstrained(P_jt[k, ]))
        }
      }
    }
  }

  theta <- c(theta_priors, theta_p_y1d, theta_Ps_pre, theta_Ps_control_post, theta_Ps_treated_post)

  structure <- list(
    J = J,
    K = K,
    T_max = T_max,
    treated_period = treated_period,
    n_pre_periods = n_pre_periods,
    n_post_periods = n_post_periods,
    len_priors = length(theta_priors),
    len_p_y1d = length(theta_p_y1d),
    len_Ps_pre = length(theta_Ps_pre),
    len_Ps_control_post = length(theta_Ps_control_post),
    len_Ps_treated_post = length(theta_Ps_treated_post)
  )

  list(theta = theta, structure = structure)
}

#' Unvectorize Unconstrained Vector Back to First-Stage Parameters
#'
#' Inverse of \code{vectorize_first_stage_unconstrained()}. Reconstructs
#' first-stage transition model parameters from an unconstrained parameter vector.
#'
#' @param theta Numeric vector of unconstrained parameters.
#' @param structure Metadata from \code{vectorize_first_stage_unconstrained()}.
#' @return A list with: \code{priors}, \code{p_y1d}, \code{Ps_control}, \code{Ps_treated}.
#' @keywords internal
#' @export
unvectorize_first_stage_unconstrained <- function(theta, structure) {
  J <- structure$J
  K <- structure$K
  n_pre_periods <- structure$n_pre_periods
  n_post_periods <- structure$n_post_periods

  idx <- 1  # Running index into theta

  # Helper to consume n elements from theta
  consume <- function(n) {
    if (n == 0) return(numeric(0))
    vals <- theta[idx:(idx + n - 1)]
    idx <<- idx + n
    vals
  }

  # 1. Priors
  priors <- unconstrained_to_simplex(consume(J - 1))

  # 2. p_y1d
  p_y1d <- lapply(1:J, function(j) {
    eta_j <- consume(2 * K - 1)
    p_flat <- unconstrained_to_simplex(eta_j)
    matrix(p_flat, nrow = K, ncol = 2)  # Column-major: first K = control, next K = treated
  })

  # 3. Ps_pre (shared pretreatment transitions)
  Ps_pre <- NULL
  if (n_pre_periods > 0) {
    Ps_pre <- lapply(1:J, function(j) {
      lapply(1:n_pre_periods, function(t) {
        P_jt <- matrix(0, nrow = K, ncol = K)
        for (k in 1:K) {
          P_jt[k, ] <- unconstrained_to_simplex(consume(K - 1))
        }
        P_jt
      })
    })
  }

  # 4. Ps_control_post
  Ps_control_post <- NULL
  if (n_post_periods > 0) {
    Ps_control_post <- lapply(1:J, function(j) {
      lapply(1:n_post_periods, function(t) {
        P_jt <- matrix(0, nrow = K, ncol = K)
        for (k in 1:K) {
          P_jt[k, ] <- unconstrained_to_simplex(consume(K - 1))
        }
        P_jt
      })
    })
  }

  # 5. Ps_treated_post
  Ps_treated_post <- NULL
  if (n_post_periods > 0) {
    Ps_treated_post <- lapply(1:J, function(j) {
      lapply(1:n_post_periods, function(t) {
        P_jt <- matrix(0, nrow = K, ncol = K)
        for (k in 1:K) {
          P_jt[k, ] <- unconstrained_to_simplex(consume(K - 1))
        }
        P_jt
      })
    })
  }

  # Assemble full Ps_control and Ps_treated
  Ps_control <- concatenate_Ps(Ps_pre, Ps_control_post)
  Ps_treated <- concatenate_Ps(Ps_pre, Ps_treated_post)

  list(
    priors = priors,
    p_y1d = p_y1d,
    Ps_control = Ps_control,
    Ps_treated = Ps_treated
  )
}

# --- Objective Function ---

#' Compute Negative Log-Likelihood for First-Stage Solver
#'
#' Objective function for the nonlinear solver. Takes an unconstrained
#' parameter vector, converts to model parameters via softmax, and computes
#' the negative observed-data log-likelihood.
#'
#' @param theta Unconstrained parameter vector.
#' @param structure Metadata from \code{vectorize_first_stage_unconstrained()}.
#' @param data_for_est Data for estimation from \code{split_data_for_estimation()}.
#' @param unit_weights Unit weights vector.
#' @return Scalar negative log-likelihood (for minimization).
#' @keywords internal
#' @export
first_stage_neg_log_likelihood <- function(theta, structure, data_for_est, unit_weights) {
  # Guard against NaN/Inf parameters from solver exploration
  if (any(!is.finite(theta))) return(.Machine$double.xmax)

  # Recover model parameters (with error handling for degenerate cases)
  result <- tryCatch({
    model_params <- unvectorize_first_stage_unconstrained(theta, structure)

    # Split p_y1d into control and treated initial distributions
    p_y1d_splitted <- split_p_y1d(model_params$p_y1d)

    is_control <- data_for_est$is_control
    is_treated <- data_for_est$is_treated

    # Log-likelihood for control units
    ll_control <- transitions_to_log_likelihood(
      data_for_est$y_indices_control,
      p_y1d_splitted$p_y1d_control,
      model_params$Ps_control,
      model_params$priors,
      unit_weights[is_control]
    )

    # Log-likelihood for treated units
    ll_treated <- transitions_to_log_likelihood(
      data_for_est$y_indices_treated,
      p_y1d_splitted$p_y1d_treated,
      model_params$Ps_treated,
      model_params$priors,
      unit_weights[is_treated]
    )

    neg_ll <- -(ll_control + ll_treated)
    if (!is.finite(neg_ll)) return(.Machine$double.xmax)
    neg_ll
  }, error = function(e) {
    .Machine$double.xmax
  })

  return(result)
}

# --- Analytic Gradient ---

#' Compute Analytic Gradient of Negative Log-Likelihood
#'
#' Computes the gradient of the negative observed-data log-likelihood with
#' respect to the unconstrained parameter vector theta. Uses the EM gradient
#' identity: the observed-data score equals the expected complete-data score
#' under the posterior distribution, combined with the softmax chain rule.
#'
#' The gradient for each softmax-parameterized simplex block has the form:
#' \eqn{dL/d\eta_m = n_m^w - p_m \cdot N^w} where \eqn{n_m^w} is the
#' weighted count in category m and \eqn{N^w} is the total weighted count.
#'
#' @param theta Unconstrained parameter vector.
#' @param structure Metadata from \code{vectorize_first_stage_unconstrained()}.
#' @param data_for_est Data for estimation from \code{split_data_for_estimation()}.
#' @param unit_weights Unit weights vector.
#' @return Numeric vector of same length as theta (gradient of negative log-likelihood).
#' @keywords internal
#' @export
first_stage_neg_log_likelihood_gradient <- function(theta, structure, data_for_est, unit_weights) {
  # Guard against NaN/Inf parameters
  if (any(!is.finite(theta))) return(rep(0, length(theta)))

  result <- tryCatch({
    # 1. Recover model parameters
    model_params <- unvectorize_first_stage_unconstrained(theta, structure)
    priors <- model_params$priors
    p_y1d <- model_params$p_y1d
    Ps_control <- model_params$Ps_control
    Ps_treated <- model_params$Ps_treated

    J <- structure$J
    K <- structure$K
    n_pre_periods <- structure$n_pre_periods
    n_post_periods <- structure$n_post_periods

    is_control <- data_for_est$is_control
    is_treated <- data_for_est$is_treated

    # 2. Split p_y1d and compute posteriors
    p_y1d_splitted <- split_p_y1d(p_y1d)

    posteriors_control <- transitions_to_posteriors(
      data_for_est$y_indices_control,
      p_y1d_splitted$p_y1d_control,
      Ps_control,
      priors
    )

    posteriors_treated <- transitions_to_posteriors(
      data_for_est$y_indices_treated,
      p_y1d_splitted$p_y1d_treated,
      Ps_treated,
      priors
    )

    # 3. Compute weighted posteriors: wZ[i,j] = unit_weight_i * posterior_ij
    w_control <- unit_weights[is_control]
    w_treated <- unit_weights[is_treated]
    wZ_control <- posteriors_control * w_control  # N_c x J
    wZ_treated <- posteriors_treated * w_treated  # N_t x J

    # Ensure matrix form even when J=1
    wZ_control <- matrix(wZ_control, ncol = J)
    wZ_treated <- matrix(wZ_treated, ncol = J)

    y_control <- data_for_est$y_indices_control
    y_treated <- data_for_est$y_indices_treated

    grad <- numeric(0)

    # --- Block A: Priors (J-1 values) ---
    if (J > 1) {
      Z_bar <- colSums(wZ_control) + colSums(wZ_treated)
      grad <- softmax_gradient(Z_bar, priors)
    }

    # --- Block B: p_y1d (J blocks of 2K-1 values) ---
    for (j in 1:J) {
      # Weighted counts for initial outcomes
      counts_control <- numeric(K)
      counts_treated <- numeric(K)
      for (k in 1:K) {
        counts_control[k] <- sum(wZ_control[y_control[, 1] == k, j])
        counts_treated[k] <- sum(wZ_treated[y_treated[, 1] == k, j])
      }
      counts_2K <- c(counts_control, counts_treated)
      p_flat <- c(p_y1d[[j]])  # Column-major: first K = control, next K = treated
      grad <- c(grad, softmax_gradient(counts_2K, p_flat))
    }

    # --- Block C: Ps_pre (shared pretreatment transitions) ---
    if (n_pre_periods > 0) {
      for (j in 1:J) {
        for (t in 1:n_pre_periods) {
          # Combine control and treated for pretreatment (shared transitions)
          y_from <- c(y_control[, t], y_treated[, t])
          y_to <- c(y_control[, t + 1], y_treated[, t + 1])
          w_j <- c(wZ_control[, j], wZ_treated[, j])
          counts <- compute_weighted_transition_counts(y_from, y_to, w_j, K)
          for (k in 1:K) {
            grad <- c(grad, softmax_gradient(counts[k, ], Ps_control[[j]][[t]][k, ]))
          }
        }
      }
    }

    # --- Block D: Ps_control_post (post-treatment control transitions) ---
    if (n_post_periods > 0) {
      for (j in 1:J) {
        for (t_idx in 1:n_post_periods) {
          t <- n_pre_periods + t_idx  # Index into full Ps list
          y_from <- y_control[, t]
          y_to <- y_control[, t + 1]
          counts <- compute_weighted_transition_counts(y_from, y_to, wZ_control[, j], K)
          for (k in 1:K) {
            grad <- c(grad, softmax_gradient(counts[k, ], Ps_control[[j]][[t]][k, ]))
          }
        }
      }
    }

    # --- Block E: Ps_treated_post (post-treatment treated transitions) ---
    if (n_post_periods > 0) {
      for (j in 1:J) {
        for (t_idx in 1:n_post_periods) {
          t <- n_pre_periods + t_idx  # Index into full Ps list
          y_from <- y_treated[, t]
          y_to <- y_treated[, t + 1]
          counts <- compute_weighted_transition_counts(y_from, y_to, wZ_treated[, j], K)
          for (k in 1:K) {
            grad <- c(grad, softmax_gradient(counts[k, ], Ps_treated[[j]][[t]][k, ]))
          }
        }
      }
    }

    stopifnot(length(grad) == length(theta))
    -grad  # Negate for minimization of negative log-likelihood

  }, error = function(e) {
    rep(0, length(theta))
  })

  return(result)
}

# --- Option Getters ---

#' Retrieve the Long Run Method for Two-Stage Multistart
#'
#' Returns the optimization method to use in the long-run stage:
#' \code{"em"} (default) for EM algorithm, or \code{"solver"} for
#' direct likelihood maximization via \code{stats::nlminb}.
#'
#' @param opts A list of options.
#' @return Character: \code{"em"} or \code{"solver"}.
#' @keywords internal
#' @export
get_first_stage_long_run_method <- function(opts) {
  get_option(opts, "long_run_method", "solver")
}

#' Retrieve Maximum Iterations for Solver
#'
#' @param opts A list of options.
#' @return Integer maximum iterations for the solver.
#' @keywords internal
get_solver_max_iter <- function(opts) {
  get_option(opts, "solver_max_iter", 500)
}

#' Retrieve Relative Tolerance for Solver
#'
#' @param opts A list of options.
#' @return Numeric relative tolerance for the solver.
#' @keywords internal
get_solver_rel_tol <- function(opts) {
  get_option(opts, "solver_rel_tol", 1e-6)
}

#' Retrieve Whether to Use Analytic Gradients in Solver
#'
#' @param opts A list of options.
#' @return Logical; \code{TRUE} (default) to use analytic gradients,
#'   \code{FALSE} to fall back to finite-difference approximation.
#' @keywords internal
#' @export
get_solver_use_gradient <- function(opts) {
  get_option(opts, "solver_use_gradient", TRUE)
}

# --- Solver Wrapper ---

#' Run Nonlinear Solver for First-Stage Estimation
#'
#' Takes an EM-converged transition model as starting point and refines
#' it using \code{stats::nlminb} with softmax reparameterization. Returns
#' output in the same format as \code{first_stage_from_guess()}.
#'
#' @param transition_model_start EM-converged transition model (starting point).
#' @param data_for_est Data structure from \code{split_data_for_estimation()}.
#' @param unit_weights Unit weights vector, or \code{NULL} for equal weights.
#' @param opts Algorithm options list.
#' @return A list with:
#'   \describe{
#'     \item{transition_model}{The optimized transition model.}
#'     \item{posteriors}{Posterior probabilities (N x J matrix).}
#'   }
#' @export
first_stage_from_solver <- function(transition_model_start, data_for_est,
                                     unit_weights = NULL, opts = list()) {
  # Validate unit_weights
  unit_weights <- check_unit_weights_validity(
    posteriors = data_for_est$y_indices_matrix,
    unit_weights = unit_weights
  )

  # Vectorize starting model to unconstrained theta
  vec_result <- vectorize_first_stage_unconstrained(transition_model_start, data_for_est)
  theta_start <- vec_result$theta
  structure <- vec_result$structure

  # Solver parameters
  solver_max_iter <- get_solver_max_iter(opts)
  solver_rel_tol <- get_solver_rel_tol(opts)
  use_gradient <- get_solver_use_gradient(opts)
  gradient_fn <- if (use_gradient) first_stage_neg_log_likelihood_gradient else NULL

  # Run nlminb
  optim_result <- stats::nlminb(
    start = theta_start,
    objective = first_stage_neg_log_likelihood,
    gradient = gradient_fn,
    structure = structure,
    data_for_est = data_for_est,
    unit_weights = unit_weights,
    control = list(
      iter.max = solver_max_iter,
      eval.max = solver_max_iter * 2,
      rel.tol = solver_rel_tol
    )
  )

  # Recover optimized model parameters
  model_params <- unvectorize_first_stage_unconstrained(optim_result$par, structure)

  # Build transition model in the same format as M_step output
  J <- structure$J
  is_control <- data_for_est$is_control
  is_treated <- data_for_est$is_treated

  transition_model <- list(
    priors = model_params$priors,
    p_y1d = model_params$p_y1d,
    Ps_control = model_params$Ps_control,
    Ps_treated = model_params$Ps_treated,
    log_likelihood = -optim_result$objective
  )

  # Compute posteriors via E-step
  posteriors <- E_step(transition_model, data_for_est)

  # Compute priors_treated from treated posteriors
  posteriors_treated <- matrix(posteriors[is_treated, ], ncol = J)
  weights_treated <- unit_weights[is_treated]
  priors_treated <- M_step_priors(posteriors_treated, weights_treated)
  transition_model$priors_treated <- priors_treated

  # Compute pmfs_initial (split from p_y1d)
  p_y1d_splitted <- split_p_y1d(model_params$p_y1d)
  transition_model$pmfs_initial_control <- p_y1d_splitted$p_y1d_control
  transition_model$pmfs_initial_treated <- p_y1d_splitted$p_y1d_treated

  # Compute empirical Ps
  Ps_from_pre <- compute_Ps_empirical(posteriors, data_for_est, unit_weights,
                                       from_0 = FALSE)
  Ps_from_0 <- compute_Ps_empirical(posteriors, data_for_est, unit_weights,
                                     from_0 = TRUE)
  transition_model <- c(transition_model, Ps_from_pre, Ps_from_0)

  return(list(transition_model = transition_model, posteriors = posteriors))
}
