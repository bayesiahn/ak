#' @title Bootstrap Functions for Transition Models
#'
#' @description
#' Functions for generating bootstrapped transition models using weighted
#' bootstrap and multiplier bootstrap methods.
#'
#' @name bootstrap
NULL

# -----------------------------------------------------------------------------
# Bootstrap Utility Functions
# -----------------------------------------------------------------------------

#' Generate Random Unit Weights
#'
#' Generates random unit weights for a specified number of units.
#' The weights are sampled from an exponential distribution and normalized
#' so that they sum to 1.
#'
#' @param N An integer specifying the number of units.
#' @return A numeric vector of length N representing the normalized unit weights.
#' @keywords internal
generate_weights <- function(N) {
  raw_weights <- rexp(N, rate = 1)
  normalized_weights <- raw_weights / sum(raw_weights)
  return(normalized_weights)
}

#' Validate and Standardize Unit Weights
#'
#' Checks if unit_weights is provided and matches the dimensions of posteriors.
#' If unit_weights is NULL, it defaults to equal weights.
#'
#' @param posteriors A matrix of posterior probabilities.
#' @param unit_weights A numeric vector of weights for each unit, or NULL.
#' @return A numeric vector of unit weights.
#' @keywords internal
check_unit_weights_validity <- function(posteriors, unit_weights) {
  if (is.null(unit_weights)) {
    unit_weights <- rep(1, nrow(posteriors))
  }
  if (length(unit_weights) != nrow(posteriors)) {
    stop("Length of unit_weights must match the number of rows in posteriors.")
  }
  return(unit_weights)
}

#' Retrieve the Number of Bootstrap Samples
#'
#' Retrieves the number of bootstrap samples to be computed based on options.
#' Default value is 200 if not specified.
#'
#' @param opts An optional list containing algorithm options.
#' @return An integer representing the number of bootstrap samples.
#' @export
get_bootstrap_counts <- function(opts) {
  get_option(opts, "bootstrap_counts", 200)
}

#' Retrieve Multistart Counts for Bootstrap Estimation
#'
#' Retrieves the number of multistart points to use during bootstrap estimation.
#' This is separate from the main estimation's multistart_counts to allow
#' fewer starting points in bootstrap (since the original estimate is used as
#' an initial guess). Default value is 10.
#'
#' @param opts An optional list containing algorithm options.
#' @return An integer representing the number of multistart points for bootstrap.
#' @export
get_multistart_bootstrap_counts <- function(opts) {
  get_option(opts, "multistart_bootstrap_counts", 10)
}

#' Retrieve Standard Error Method Option
#'
#' Retrieves the standard error method specified in the options.
#' Returns "bootstrap", "bootstrap_multiplier", or "none".
#'
#' @param opts An optional list containing algorithm options.
#' @return A character string indicating the SE method.
#' @export
get_se_method <- function(opts) {
  if (is.list(opts) && "se_method" %in% names(opts)) {
    if (opts$se_method == "boot" || opts$se_method == "bootstrap") {
      return("bootstrap")
    }
    if (opts$se_method == "multiplier" ||
        opts$se_method == "boot_multiplier" ||
        opts$se_method == "bootstrap_multiplier") {
      return("bootstrap_multiplier")
    }
  }
  return("none")
}

#' Retrieve the Number of Cores for Parallel Bootstrap
#'
#' Retrieves the number of cores to be used in parallel computation.
#' Defaults to min(10, detectCores() - 1) if not specified.
#'
#' @param opts An optional list containing algorithm options.
#' @return An integer representing the number of cores.
#' @export
get_bootstrap_mc_core_counts <- function(opts) {
  default_cores <- min(10, max(1, parallel::detectCores() - 1))
  get_option(opts, "bootstrap_mc_core_counts", default_cores)
}

#' Generate Rademacher Distributed Multipliers
#'
#' Generates a vector of N random multipliers from the Rademacher distribution,
#' where each value is either -1 or 1 with equal probability.
#'
#' @param N An integer representing the number of multipliers to generate.
#' @return A numeric vector of length N with values -1 or 1.
#' @keywords internal
draw_multipliers <- function(N) {
  sample(c(-1, 1), N, replace = TRUE)
}

# -----------------------------------------------------------------------------
# Bootstrap Generation Functions
# -----------------------------------------------------------------------------

#' Generate Bootstrapped Transition Models from a Wide Matrix
#'
#' Generates bootstrapped transition models by repeatedly resampling the input
#' data and applying the estimate_transition_model_from_wide_matrix function
#' to each resample. Uses parallel processing with mclapply.
#'
#' @param J The number of latent groups in the model.
#' @param y A matrix or data frame of outcomes for each unit over time.
#' @param g A vector indicating the treatment status for each unit.
#' @param c (Optional) A vector indicating the cluster for each unit.
#' @param transition_model_guesses (Optional) A list of initial guesses.
#' @param outcomes_of_interest (Optional) A vector of outcome names.
#' @param opts A list of options including:
#'   \itemize{
#'     \item bootstrap_counts: Number of bootstrap resamples (default 200)
#'     \item mc_cores: Number of cores for parallel processing (default 4)
#'   }
#'
#' @return A list of bootstrapped transition models.
#' @export
bootstrap_transition_models_from_wide_matrix <- function(J, y, g, c = NULL,
                                                         transition_model_guesses = NULL,
                                                         outcomes_of_interest = NULL,
                                                         opts = NULL) {
  bootstrap_counts <- get_bootstrap_counts(opts)

  # Assign cluster identifiers
  if (is.null(c)) {
    c <- seq_along(g)
  }
  n_clusters <- length(unique(c))
  c_factor <- as.factor(c)

  # Generate seeds for reproducibility (one per bootstrap sample)
  # This avoids storing all weights upfront - ~99% memory reduction
  set.seed(get_option(opts, "bootstrap_seed", 1234))
  bootstrap_seeds <- sample.int(.Machine$integer.max, bootstrap_counts)

  # Prepare bootstrap options: disable multistart parallelization inside
  # bootstrap to avoid nested fork issues
  opts_bootstrap <- opts
  if (is.null(opts_bootstrap)) opts_bootstrap <- list()
  opts_bootstrap$parallel_multistart <- FALSE
  # Use multistart_bootstrap_counts instead of multistart_counts for bootstrap
  # This allows fewer starting points since the original estimate is used as an initial guess
  opts_bootstrap$multistart_counts <- get_multistart_bootstrap_counts(opts)
  # Set two-stage multistart parameters for bootstrap
  # Use same value for short run counts (explicit)
  opts_bootstrap$multistart_short_run_counts <- get_multistart_bootstrap_counts(opts)
  # Use fewer long run candidates since we start from original estimate
  opts_bootstrap$multistart_long_run_counts <- 4

  # Estimate transition models for each bootstrap sample
  # Weights are generated on-the-fly to save memory
  mc_core_counts <- get_bootstrap_mc_core_counts(opts)
  bootstrap_estimates <- parallel::mclapply(1:bootstrap_counts, function(b) {
    # Generate weights on-the-fly using pre-generated seed
    set.seed(bootstrap_seeds[b])
    cluster_weights <- generate_weights(n_clusters)
    unit_weights <- cluster_weights[c_factor]

    estimate_transition_model_from_wide_matrix(
      J = J,
      y = y,
      g = g,
      outcomes_of_interest = outcomes_of_interest,
      unit_weights = unit_weights,
      transition_model_guesses = transition_model_guesses,
      opts = opts_bootstrap
    )
  }, mc.cores = mc_core_counts)

  # Filter out NULL results from failed mclapply forks (e.g., segfaults)
  n_total <- length(bootstrap_estimates)
  bootstrap_estimates <- Filter(Negate(is.null), bootstrap_estimates)
  n_failed <- n_total - length(bootstrap_estimates)
  if (n_failed > 0) {
    warning(sprintf("%d of %d bootstrap replications failed and were dropped.", n_failed, n_total))
  }

  return(bootstrap_estimates)
}

#' Perform Multiplier Bootstrap for Model Estimation
#'
#' Performs a multiplier bootstrap to estimate a transition model by resampling
#' with Rademacher multipliers. Bootstrap estimates are generated by applying
#' randomly drawn multipliers to model scores and reconstructing the model.
#'
#' @param y A matrix of observed outcomes (rows = individuals, cols = time).
#' @param g A numeric vector indicating treatment assignment.
#' @param model A list containing the transition model components.
#' @param opts An optional list of bootstrap options with bootstrap_counts.
#'
#' @return A list with:
#'   \itemize{
#'     \item bootstrap_estimates: List of estimated models from each iteration
#'     \item scores: The computed scores from the original model
#'     \item bootstrap_multipliers: List of multipliers used
#'   }
#' @export
multiplier_bootstrap <- function(y, g, model, opts = NULL) {
  bootstrap_counts <- get_bootstrap_counts(opts)
  N <- length(g)

  # Compute scores from the model
  scores <- model_to_scores(y, g, model)

  # Draw multipliers and corresponding models
  bootstrap_estimates <- list()
  multipliers <- list()
  for (i in 1:bootstrap_counts) {
    multipliers_i <- draw_multipliers(N)
    model_bootstrap <- scores_to_model(scores, multipliers_i)
    bootstrap_estimates[[i]] <- model_bootstrap
    multipliers[[i]] <- multipliers_i
  }

  list(bootstrap_estimates = bootstrap_estimates, scores = scores,
       bootstrap_multipliers = multipliers)
}
