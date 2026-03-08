#' Retrieve Maximum Iteration Count for the First Stage
#'
#' This function retrieves the maximum number of iterations to be used in the
#' first stage of an algorithm. If not specified in the options, a default
#' value is returned.
#'
#' @param opts An optional list containing the algorithm's options, specifically
#'        the maximum number of iterations under `$max_iter`.
#'
#' @return The maximum number of iterations to use.
#'
#' @examples
#' \dontrun{
#' opts <- list(max_iter = 200)
#' max_iter <- get_first_stage_max_iter(opts)
#' # Default value
#' max_iter_default <- get_first_stage_max_iter(NULL)
#' }
get_first_stage_max_iter <- function(opts) {
  get_option(opts, "max_iter", 100)
}

#' Retrieve Tolerance Level for the First Stage
#'
#' This function gets the tolerance level for convergence in the first stage of
#' an algorithm. If not specified in the options, a default value is returned.
#'
#' @param opts An optional list containing the algorithm's options, specifically
#'        the tolerance level under `$tol`.
#'
#' @return The tolerance level for convergence.
#'
#' @examples
#' \dontrun{
#' opts <- list(tol = 1e-5)
#' tol <- get_first_stage_tol(opts)
#' # Default value
#' tol_default <- get_first_stage_tol(NULL)
#' }
get_first_stage_tol <- function(opts) {
  get_option(opts, "tol", 1e-3)
}

#' Retrieve Multistart Counts for the First Stage
#'
#' This function gets the number of multistarts to use in the first stage of an
#' algorithm. If not specified in the options, a default value proportional to
#' the number of latent groups (J) is returned.
#'
#' @param opts An optional list containing the algorithm's options, specifically
#'       the number of multistarts under `$multistart_counts` or a base count
#'       under `$multistart_counts_base`.
#' @param J The number of latent groups. The default multistart count scales
#'       linearly with J (base_count * J).
#'
#' @return The number of multistarts to use.
#'
#' @examples
#' \dontrun{
#' opts <- list(multistart_counts = 20)
#' multistart_counts <- get_first_stage_multistart_counts(opts, J = 2)
#' # Default value scales with J: J=2 -> 20, J=3 -> 40, ..., capped at 100
#' multistart_counts_default <- get_first_stage_multistart_counts(NULL, J = 3)  # Returns 40
#' }
get_first_stage_multistart_counts <- function(opts, J = 1) {
  # If explicit multistart_counts is provided, use it directly
  explicit_count <- get_option(opts, "multistart_counts", NULL)
  if (!is.null(explicit_count)) {
    return(explicit_count)
  }

  # For J=1, no latent group structure to estimate, so only 1 start needed
  if (J == 1) {
    return(1)
  }

  # Otherwise, scale proportionally with (J-1): J=2 -> 20, J=3 -> 40, etc.
  # Capped at max_count (default 100) to prevent excessive computation
  base_count <- get_option(opts, "multistart_counts_base", 2000)

  max_count <- get_option(opts, "multistart_counts_max", 4000)
  min(base_count * (J - 1), max_count)
}

#' Retrieve Short Run Multistart Counts
#'
#' Returns the number of initial guesses to use in the short run stage of
#' two-stage multistart optimization. Defaults to the standard multistart_counts
#' if not specified.
#'
#' @param opts A list of options, may contain `multistart_short_run_counts`.
#' @param J The number of latent groups.
#' @return Integer number of short run guesses.
#' @keywords internal
get_first_stage_short_run_counts <- function(opts, J = 1) {
  explicit_count <- get_option(opts, "multistart_short_run_counts", NULL)
  if (!is.null(explicit_count)) {
    return(explicit_count)
  }
  # Default to standard multistart_counts
  get_first_stage_multistart_counts(opts, J)
}

#' Retrieve Long Run Multistart Counts
#'
#' Returns the number of top candidates to keep for the long run stage of
#' two-stage multistart optimization. Defaults to 4 if not specified.
#'
#' @param opts A list of options, may contain `multistart_long_run_counts`.
#' @return Integer number of candidates to keep for long run.
#' @keywords internal
get_first_stage_long_run_counts <- function(opts) {
  get_option(opts, "multistart_long_run_counts", 20)
}

#' Retrieve Maximum Iterations for Short Run Stage
#'
#' Returns the maximum number of EM iterations for the short run stage.
#' Defaults to 20 if not specified.
#'
#' @param opts A list of options, may contain `max_iter_short_run`.
#' @return Integer maximum iterations for short run.
#' @keywords internal
get_first_stage_max_iter_short_run <- function(opts) {
  get_option(opts, "max_iter_short_run", 200)
}

#' Retrieve Prior Lower Bound
#'
#' Returns the minimum acceptable prior probability for any latent group.
#' Models with any prior below this threshold are filtered out during
#' multistart selection to prevent degenerate mixture components.
#' Default is 0 (no filtering).
#'
#' @param opts A list of options, may contain \code{prior_lower_bound}.
#' @return Numeric lower bound for priors.
#' @keywords internal
get_prior_lower_bound <- function(opts) {
  get_option(opts, "prior_lower_bound", 0)
}

#' Retrieve Maximum Iterations for Long Run Stage
#'
#' Returns the maximum number of EM iterations for the long run stage.
#' Defaults to 200 if not specified.
#'
#' @param opts A list of options, may contain `max_iter_long_run`.
#' @return Integer maximum iterations for long run.
#' @keywords internal
get_first_stage_max_iter_long_run <- function(opts) {
  get_option(opts, "max_iter_long_run", 200)
}

#' Check if Two-Stage Multistart is Enabled
#'
#' Determines whether to use the two-stage multistart optimization approach.
#' This is enabled by default when multistart_counts > 1, but can be disabled
#' by setting `two_stage_multistart = FALSE` in opts.
#'
#' @param opts A list of options, may contain `two_stage_multistart`.
#' @return Logical indicating whether to use two-stage multistart.
#' @keywords internal
use_two_stage_multistart <- function(opts) {
  get_option(opts, "two_stage_multistart", TRUE)
}

#' Check if Multistart Should Run in Parallel
#'
#' Determines whether multistart estimation should use parallel processing.
#' Parallelization is disabled when: (1) explicitly set to FALSE in opts,
#' (2) running inside another parallel context (e.g., bootstrap), or
#' (3) only one core is available.
#'
#' @param opts A list of options, may contain `parallel_multistart`.
#' @return Logical indicating whether to use parallel multistart.
#' @keywords internal
should_parallelize_multistart <- function(opts) {
 # Check explicit option
 parallel_opt <- get_option(opts, "parallel_multistart", "auto")

 # If explicitly disabled, return FALSE
 if (identical(parallel_opt, FALSE)) {
   return(FALSE)
 }

 # If explicitly enabled, return TRUE
 if (identical(parallel_opt, TRUE)) {
   return(TRUE)
 }

 # Auto mode: check available cores
 available_cores <- get_multistart_mc_core_counts(opts)
 if (available_cores <= 1) {
   return(FALSE)
 }

 return(TRUE)
}

#' Retrieve the Number of Cores for Parallel Multistart
#'
#' Retrieves the number of cores to use for parallel multistart.
#' Defaults to min(10, detectCores() - 1) if not specified.
#'
#' @param opts An optional list containing algorithm options.
#' @return An integer representing the number of cores.
#' @keywords internal
get_multistart_mc_core_counts <- function(opts) {
 default_cores <- min(10, max(1, parallel::detectCores() - 1))
 get_option(opts, "multistart_mc_core_counts", default_cores)
}

#' Retrieve the Number of Cores for Short Run Stage
#'
#' Retrieves the number of cores to use for the short run stage of two-stage multistart.
#' Defaults to `multistart_mc_core_counts` if `multistart_mc_core_counts_short_run` is not specified.
#'
#' @param opts An optional list containing algorithm options.
#' @return An integer representing the number of cores for short run.
#' @keywords internal
get_multistart_mc_core_counts_short_run <- function(opts) {
  # First check for short-run specific setting
  short_run_cores <- get_option(opts, "multistart_mc_core_counts_short_run", NULL)
  if (!is.null(short_run_cores)) {
    return(short_run_cores)
  }
  # Fall back to general multistart_mc_core_counts
  get_multistart_mc_core_counts(opts)
}

#' Retrieve the Number of Cores for Long Run Stage
#'
#' Retrieves the number of cores to use for the long run stage of two-stage multistart.
#' Defaults to `multistart_mc_core_counts` if `multistart_mc_core_counts_long_run` is not specified.
#'
#' @param opts An optional list containing algorithm options.
#' @return An integer representing the number of cores for long run.
#' @keywords internal
get_multistart_mc_core_counts_long_run <- function(opts) {
  # First check for long-run specific setting
  long_run_cores <- get_option(opts, "multistart_mc_core_counts_long_run", NULL)
  if (!is.null(long_run_cores)) {
    return(long_run_cores)
  }
  # Fall back to general multistart_mc_core_counts
  get_multistart_mc_core_counts(opts)
}

#' Run the First Stage of an Algorithm from an Initial Guess
#'
#' This function executes the first stage of an algorithm, starting from a
#' given initial transition model guess. It iterates through E and M steps until
#' convergence or the maximum number of iterations is reached.
#'
#' @param transition_model_guess An initial guess for the transition model.
#' @param data_for_est A list returned by \code{split_data_for_estimation()}, containing
#'        state sequences and indicators split by treatment status.
#' @param opts (Optional) A list of options for the algorithm, including
#'        maximum iterations and tolerance level.
#' @param unit_weights A numeric vector of weights for each unit, or \code{NULL}.
#'        If \code{NULL}, equal weights are assumed.
#' @param max_iter_override (Optional) If provided, overrides the max_iter from opts.
#'        Used internally by two-stage multistart to control short/long run iterations.
#' @param compute_empirical_Ps Logical; if TRUE (default), computes empirical transition
#'        matrices from pre-treatment periods. Set to FALSE during short run to save computation.
#'
#' @return The transition model after running the first stage of the algorithm.
#'
#' @export
first_stage_from_guess <- function(transition_model_guess, data_for_est,
                                   opts = list(),
                                   unit_weights = NULL,
                                   max_iter_override = NULL,
                                   compute_empirical_Ps = TRUE) {
  # Check inputs
  max_iter <- if (!is.null(max_iter_override)) max_iter_override else get_first_stage_max_iter(opts)
  tol <- get_first_stage_tol(opts)
  unit_weights <- check_unit_weights_validity(posteriors = data_for_est$y_indices_matrix, unit_weights)

  # Initialize transition model
  transition_model <- transition_model_guess

  # Initialize log likelihood
  log_likelihood <- -Inf

  # Iterate until convergence
  for (iter in 1:max_iter) {
    # E-Step
    posteriors <- E_step(transition_model, data_for_est)

    # M-Step
    transition_model <- M_step(posteriors, data_for_est,
                         unit_weights = unit_weights)

    # Check convergence using relative change for numerical stability
    new_log_likelihood <- transition_model$log_likelihood
    rel_change <- abs(new_log_likelihood - log_likelihood) /
                  (1 + abs(new_log_likelihood))
    if (rel_change < tol) break
    log_likelihood <- new_log_likelihood
  }

  # Compute empirical Ps only if requested (skip during short run for efficiency)
  if (compute_empirical_Ps) {
    Ps_from_pre <- compute_Ps_empirical(posteriors, data_for_est, unit_weights,
                                        from_0 = FALSE)
    Ps_from_0 <- compute_Ps_empirical(posteriors, data_for_est, unit_weights,
                                      from_0 = TRUE)
    transition_model <- c(transition_model, Ps_from_pre, Ps_from_0)
  }

  # Return transition model with posteriors
  return(list(transition_model = transition_model, posteriors = posteriors))
}

#' Generate an Initial Guess for the First Stage of an Algorithm
#'
#' This function generates an initial guess for the first stage of an algorithm
#' based on the provided state sequences (`y_indices_matrix`) and the number
#' of latent groups (`J`). It creates a transition model with specific outcomes and
#' time periods.
#'
#' @param y_indices_matrix A matrix of state sequences for each unit.
#' @param J The number of latent groups in the model.
#' @param treated_period (Optional) The period when treatment starts. Defaults to 0,
#'
#' @return An initial transition model guess, with specified outcomes and number of latent groups.
#'
#' @examples
#' \dontrun{
#' y_indices_matrix <- matrix(sample(1:3, 9, replace = TRUE), nrow = 3)
#' J <- 2
#' initial_guess <- generate_first_stage_guess(y_indices_matrix, J)
#' }
generate_first_stage_guess <- function(y_indices_matrix, J, treated_period = 0) {
  # Extract unique outcomes from the y_indices_matrix
  outcomes <- sort(unique(c(y_indices_matrix)))

  # Generate the initial transition model based on the outcomes, number of latent groups (J),
  # and the number of time periods (T_max) derived from the y_indices_matrix
  # NOTE: treated_period is set to be zero because nobody is treated in pre-treatment periods
  generate_transition_model(outcomes = outcomes, J = J,
                      T_max = ncol(y_indices_matrix),
                      treated_period = treated_period)
}

#' Select Top Candidates by Log-Likelihood
#'
#' Selects the top K models by log-likelihood from a list of estimated models.
#'
#' @param estimated_models A list of estimated models, each with a `transition_model`
#'        element containing a `log_likelihood` field.
#' @param k The number of top candidates to select.
#' @return A list of the top K models.
#' @keywords internal
select_top_candidates <- function(estimated_models, k, opts = list()) {
  # Filter out failed results (NULL, errors, or non-list objects)
  valid_mask <- sapply(estimated_models, function(x) {
    is.list(x) && !is.null(x$transition_model) &&
      !is.null(x$transition_model$log_likelihood)
  })
  valid_models <- estimated_models[valid_mask]

  if (length(valid_models) == 0) {
    stop("All multistart attempts failed. Check data or model specification.")
  }

  # Filter corner solutions: models with boundary priors or degenerate transitions
  boundary_tol <- get_option(opts, "corner_solution_tol", 0.025)
  corner_free <- sapply(valid_models, function(x) {
    tm <- x$transition_model
    # Check priors: no group should have near-zero or near-one probability
    priors <- tm$priors
    if (!is.null(priors) && length(priors) > 1) {
      if (any(priors < boundary_tol) || any(priors > 1 - boundary_tol)) {
        return(FALSE)
      }
    }
    # Check transition matrices: no row should be degenerate (all mass on one state)
    Ps <- tm$Ps_control
    if (!is.null(Ps)) {
      for (j in seq_along(Ps)) {
        for (t in seq_along(Ps[[j]])) {
          P <- Ps[[j]][[t]]
          if (any(P > 1 - boundary_tol & P < 1) || any(P < boundary_tol & P > 0)) {
            # Has near-boundary entries; check if any row is degenerate
            for (row_idx in seq_len(nrow(P))) {
              if (max(P[row_idx, ]) > 1 - boundary_tol) {
                return(FALSE)
              }
            }
          }
        }
      }
    }
    TRUE
  })

  if (sum(corner_free) > 0) {
    valid_models <- valid_models[corner_free]
  }
  # If ALL have corner solutions, keep all (fallback)

  likelihoods <- sapply(valid_models, function(x) x$transition_model$log_likelihood)
  # Handle case where k is larger than the number of valid models
  k <- min(k, length(valid_models))
  top_indices <- order(likelihoods, decreasing = TRUE)[1:k]
  valid_models[top_indices]
}

#' Perform the First Stage of an Algorithm with Multiple Initial Guesses
#'
#' This function executes the first stage of an algorithm using multiple initial
#' guesses. By default, it uses a two-stage optimization approach:
#' \enumerate{
#'   \item \strong{Short run stage:} Run EM for a small number of iterations
#'         (default 20) on many initial guesses to quickly filter out poor candidates.
#'   \item \strong{Long run stage:} Select top candidates by likelihood and run EM
#'         for more iterations (default 500) to convergence.
#' }
#'
#' The two-stage approach can be disabled by setting \code{two_stage_multistart = FALSE}
#' in opts, which reverts to the original single-stage behavior.
#'
#' @param data_for_est A list returned by \code{split_data_for_estimation()}, containing
#'        state sequences and indicators split by treatment status.
#'        It should include `y_indices_pretreatment` for the pre-treatment data.
#' @param g N-length vector indicating treatment status; ith element > 0 if ith unit is treated.
#' @param J The number of latent groups in the model.
#' @param unit_weights A numeric vector of weights for each unit, or \code{NULL}.
#'        If \code{NULL}, equal weights are assumed.
#' @param transition_model_guesses (Optional) A list of initial transition model guesses.
#'        If \code{NULL}, initial guesses are generated internally.
#' @param opts (Optional) A list of options for the algorithm. Relevant options include:
#'   \describe{
#'     \item{\code{two_stage_multistart}}{Logical; if TRUE (default), uses two-stage
#'           optimization. Set to FALSE for single-stage behavior.}
#'     \item{\code{multistart_counts}}{Number of initial guesses for single-stage mode,
#'           or legacy option. Default scales with J.}
#'     \item{\code{multistart_short_run_counts}}{Number of initial guesses for short run
#'           stage. Defaults to \code{multistart_counts}.}
#'     \item{\code{multistart_long_run_counts}}{Number of top candidates to keep for
#'           long run stage. Default is 20.}
#'     \item{\code{max_iter_short_run}}{Maximum EM iterations for short run stage.
#'           Default is 20.}
#'     \item{\code{max_iter_long_run}}{Maximum EM iterations for long run stage.
#'           Default is 500.}
#'     \item{\code{max_iter}}{Maximum EM iterations for single-stage mode. Default is 100.}
#'   }
#'
#' @return The transition model with the highest log likelihood among the candidates.
#' @export
first_stage <- function(data_for_est, g, J, unit_weights = NULL, transition_model_guesses = NULL,
                        opts = list()) {
  y_indices <- data_for_est$y_indices_matrix
  treated_period <- data_for_est$treated_period

  # Check if y_indices_pretreatment has at least two columns
  if (ncol(y_indices) < 2) {
    stop("y_indices must have at least two columns.")
  }

  # Determine whether to use two-stage multistart
  use_two_stage <- use_two_stage_multistart(opts) && J > 1

  if (use_two_stage) {
    return(first_stage_two_stage(data_for_est, g, J, unit_weights,
                                  transition_model_guesses, opts))
  } else {
    return(first_stage_single_stage(data_for_est, g, J, unit_weights,
                                     transition_model_guesses, opts))
  }
}

#' Single-Stage First Stage Estimation (Original Behavior)
#'
#' Runs the EM algorithm for max_iter iterations on all initial guesses and
#' selects the best model. This is the original behavior before two-stage
#' optimization was introduced.
#'
#' @inheritParams first_stage
#' @return The transition model with the highest log likelihood.
#' @keywords internal
first_stage_single_stage <- function(data_for_est, g, J, unit_weights = NULL,
                                      transition_model_guesses = NULL, opts = list()) {
  y_indices <- data_for_est$y_indices_matrix
  treated_period <- data_for_est$treated_period

  # Generate initial guesses if not provided, or append additional guesses
  # to reach the target multistart_counts
  multistart_counts <- get_first_stage_multistart_counts(opts, J)
  if (is.null(transition_model_guesses)) {
    transition_model_guesses <- list()
  }
  n_provided <- length(transition_model_guesses)
  if (n_provided < multistart_counts) {
    set.seed(1234)
    additional_guesses <- lapply(seq_len(multistart_counts - n_provided), function(x) {
      generate_first_stage_guess(y_indices, J, treated_period = treated_period)
    })
    transition_model_guesses <- c(transition_model_guesses, additional_guesses)
  }

  # Compute first stage for each guess with unit weights
  # Use parallel processing if enabled and not inside another parallel context
  run_first_stage <- function(x) {
    first_stage_from_guess(x, data_for_est = data_for_est, opts = opts,
                           unit_weights = unit_weights)
  }

  if (should_parallelize_multistart(opts)) {
    mc_cores <- get_multistart_mc_core_counts(opts)
    transition_model_guesses <- parallel::mclapply(
      transition_model_guesses, run_first_stage, mc.cores = mc_cores
    )
  } else {
    transition_model_guesses <- lapply(transition_model_guesses, run_first_stage)
  }

  # Filter out failed results
  valid_mask <- sapply(transition_model_guesses, function(x) {
    is.list(x) && !is.null(x$transition_model) &&
      !is.null(x$transition_model$log_likelihood)
  })
  valid_models <- transition_model_guesses[valid_mask]

  if (length(valid_models) == 0) {
    stop("All multistart attempts failed. Check data or model specification.")
  }

  # Filter degenerate priors
  prior_lb <- get_prior_lower_bound(opts)
  if (prior_lb > 0) {
    prior_ok <- sapply(valid_models, function(x) {
      priors <- x$transition_model$priors
      !is.null(priors) && all(priors >= prior_lb)
    })
    if (sum(prior_ok) > 0) {
      valid_models <- valid_models[prior_ok]
    }
  }

  # Pick one with the largest likelihood
  likelihoods <- sapply(valid_models,
                        function(x) x$transition_model$log_likelihood)
  valid_models[[which.max(likelihoods)]]
}

#' Two-Stage First Stage Estimation
#'
#' Implements the two-stage multistart optimization:
#' \enumerate{
#'   \item \strong{Short run:} Run EM for few iterations on many guesses
#'   \item \strong{Long run:} Select top candidates and run EM to convergence
#' }
#'
#' @inheritParams first_stage
#' @return The transition model with the highest log likelihood.
#' @keywords internal
first_stage_two_stage <- function(data_for_est, g, J, unit_weights = NULL,
                                   transition_model_guesses = NULL, opts = list()) {
  y_indices <- data_for_est$y_indices_matrix
  treated_period <- data_for_est$treated_period

  # Get two-stage parameters
  short_run_counts <- get_first_stage_short_run_counts(opts, J)
  long_run_counts <- get_first_stage_long_run_counts(opts)
  max_iter_short <- get_first_stage_max_iter_short_run(opts)
  max_iter_long <- get_first_stage_max_iter_long_run(opts)

  # Generate initial guesses if not provided
  if (is.null(transition_model_guesses)) {
    transition_model_guesses <- list()
  }
  n_provided <- length(transition_model_guesses)
  if (n_provided < short_run_counts) {
    set.seed(1234)
    additional_guesses <- lapply(seq_len(short_run_counts - n_provided), function(x) {
      generate_first_stage_guess(y_indices, J, treated_period = treated_period)
    })
    transition_model_guesses <- c(transition_model_guesses, additional_guesses)
  }

  # --- Stage 1: Short run ---
  # Run EM for a small number of iterations to filter out poor candidates
  run_short <- function(x) {
    first_stage_from_guess(x, data_for_est = data_for_est, opts = opts,
                           unit_weights = unit_weights,
                           max_iter_override = max_iter_short,
                           compute_empirical_Ps = FALSE)
  }

  if (should_parallelize_multistart(opts)) {
    mc_cores <- get_multistart_mc_core_counts_short_run(opts)
    short_run_results <- parallel::mclapply(
      transition_model_guesses, run_short, mc.cores = mc_cores
    )
  } else {
    short_run_results <- lapply(transition_model_guesses, run_short)
  }

  # --- Select top candidates ---
  top_candidates <- select_top_candidates(short_run_results, long_run_counts, opts)

  # --- Stage 2: Long run ---
  # Run to convergence on the top candidates using either EM or nonlinear solver
  # Use the transition_model from short run as the starting point
  long_run_method <- get_first_stage_long_run_method(opts)

  if (long_run_method == "solver") {
    run_long <- function(x) {
      first_stage_from_solver(x$transition_model, data_for_est = data_for_est,
                               unit_weights = unit_weights, opts = opts)
    }
  } else {
    run_long <- function(x) {
      first_stage_from_guess(x$transition_model, data_for_est = data_for_est, opts = opts,
                             unit_weights = unit_weights,
                             max_iter_override = max_iter_long,
                             compute_empirical_Ps = TRUE)
    }
  }

  if (should_parallelize_multistart(opts)) {
    mc_cores <- get_multistart_mc_core_counts_long_run(opts)
    long_run_results <- parallel::mclapply(
      top_candidates, run_long, mc.cores = mc_cores
    )
  } else {
    long_run_results <- lapply(top_candidates, run_long)
  }

  # Filter out failed results
  valid_mask <- sapply(long_run_results, function(x) {
    is.list(x) && !is.null(x$transition_model) &&
      !is.null(x$transition_model$log_likelihood)
  })
  valid_models <- long_run_results[valid_mask]

  if (length(valid_models) == 0) {
    stop("All multistart attempts failed. Check data or model specification.")
  }

  # Filter degenerate priors (may have drifted during long run)
  prior_lb <- get_prior_lower_bound(opts)
  if (prior_lb > 0) {
    prior_ok <- sapply(valid_models, function(x) {
      priors <- x$transition_model$priors
      !is.null(priors) && all(priors >= prior_lb)
    })
    if (sum(prior_ok) > 0) {
      valid_models <- valid_models[prior_ok]
    }
  }

  # Pick the best model
  likelihoods <- sapply(valid_models,
                        function(x) x$transition_model$log_likelihood)
  valid_models[[which.max(likelihoods)]]
}

#' Compute Transition Matrices from Pre-treatment Periods
#'
#' This function estimates post-treatment transition matrices for control and
#' treated units using only pre-treatment transition probabilities as a baseline.
#' It is typically used to assess how transition patterns from the pre-treatment
#' period extend into post-treatment periods under a baseline scenario.
#'
#' @param posteriors A matrix of posterior probabilities (n × J), where each row
#'        corresponds to a unit and each column to a latent group.
#' @param data_for_est A list returned by \code{split_data_for_estimation()},
#'        containing state sequences, transition indicators, and treatment period
#'        information. Must include:
#'        \itemize{
#'          \item \code{treated_period} — integer, period when treatment starts.
#'          \item \code{is_control}, \code{is_treated} — logical vectors indicating
#'                control and treated units.
#'          \item \code{transition_indicators_control_post_baseline} — post-treatment
#'                transition indicators for control units under baseline.
#'          \item \code{transition_indicators_treated_post_baseline} — post-treatment
#'                transition indicators for treated units under baseline.
#'        }
#' @param unit_weights Optional numeric vector of weights for each unit. If \code{NULL},
#'        equal weights are assumed.
#' @param from_0 Logical; if \code{TRUE}, estimates transition matrices starting from period 1.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{Ps_control_from_pre}}{List of estimated transition matrices for control units,
#'         using pre-treatment matrices concatenated with post-treatment matrices under baseline.}
#'   \item{\code{Ps_treated_from_pre}}{List of estimated transition matrices for treated units,
#'         using pre-treatment matrices concatenated with post-treatment matrices under baseline.}
#' }
compute_Ps_empirical <- function(posteriors, data_for_est,
                                unit_weights = NULL,
                                from_0 = FALSE) {

  # Estimate transition matrices for control and treated units in post-treatment periods
  unit_weights <- check_unit_weights_validity(posteriors, unit_weights)

  # Split posteriors and weights
  J <- ncol(posteriors)
  treated_period <- data_for_est$treated_period
  is_control <- data_for_est$is_control
  is_treated <- data_for_est$is_treated
  posteriors_control <- matrix(posteriors[is_control, ], ncol = J)
  posteriors_treated <- matrix(posteriors[is_treated, ], ncol = J)
  weights_control <- unit_weights[is_control]
  weights_treated <- unit_weights[is_treated]

  if (from_0) {
    Ps_control_from_0 <- M_step_Ps(posteriors_control,
                                  data_for_est$transition_indicators_control_baseline,
                                  weights_control)
    Ps_treated_from_0 <- M_step_Ps(posteriors_treated,
                                  data_for_est$transition_indicators_treated_baseline,
                                  weights_treated)
    return(list(Ps_control_from_0 = Ps_control_from_0,
                Ps_treated_from_0 = Ps_treated_from_0))
  }

  Ps_pre <- NULL
  if (treated_period > 2) {
    Ps_pre <- lapply(1:J, function(j) lapply(1:(treated_period-2), function(x) NA))
  }

  Ps_control_from_pre <- concatenate_Ps(Ps_pre,
                                        M_step_Ps(posteriors_control,
                                                  data_for_est$transition_indicators_control_post_baseline,
                                                  weights_control))
  Ps_treated_from_pre <- concatenate_Ps(Ps_pre,
                                        M_step_Ps(posteriors_treated,
                                                  data_for_est$transition_indicators_treated_post_baseline,
                                                  weights_treated))
  list(Ps_control_from_pre = Ps_control_from_pre,
       Ps_treated_from_pre = Ps_treated_from_pre)
}
