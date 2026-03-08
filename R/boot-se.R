#' @title Bootstrap Standard Error Functions
#'
#' @description
#' Functions for computing standard errors and quantiles from bootstrapped
#' transition models.
#'
#' @name bootstrap-se
NULL

# -----------------------------------------------------------------------------
# Elementwise Computation Utilities
# -----------------------------------------------------------------------------

#' Recursive Elementwise Computation
#'
#' Generic helper function that recursively applies a function to scalar,
#' vector, matrix, or nested list parameters.
#'
#' @param values A list of values to process (scalars, vectors, matrices,
#'   or nested lists).
#' @param func The function to apply to the elements.
#' @param ... Additional arguments to pass to func.
#'
#' @return A scalar, vector, matrix, or nested list with the result.
#' @keywords internal
recursive_elementwise <- function(values, func, ...) {
  if (length(values) == 0) {
    return(values)
  }

  first_value <- values[[1]]
  if (is.list(first_value)) {
    # Recursive case: handle lists (or nested lists)
    lapply(seq_along(first_value), function(idx) {
      sub_values <- lapply(values, `[[`, idx)
      recursive_elementwise(sub_values, func, ...)
    })
  } else if (is.matrix(first_value)) {
    # Base case: handle matrices
    array_values <- array(unlist(values), dim = c(dim(first_value), length(values)))
    apply(array_values, 1:2, func, ...)
  } else if (is.vector(first_value)) {
    # Base case: handle vectors
    matrix_values <- do.call(rbind, values)
    apply(matrix_values, 2, func, ...)
  } else {
    # Base case: handle scalars
    func(unlist(values), ...)
  }
}

#' Compute Standard Errors for Nested Structures
#'
#' Recursively computes standard errors for scalar, vector, matrix, or
#' nested list parameters.
#'
#' @param values A list of values for which standard errors are computed.
#' @return A scalar, vector, matrix, or nested list of standard errors.
#' @export
elementwise_sd <- function(values) {
  recursive_elementwise(values, sd)
}

#' Compute Elementwise Quantiles for Nested Structures
#'
#' Computes quantiles element-wise for scalar, vector, matrix, or nested
#' list parameters.
#'
#' @param values A list of values for which quantiles are computed.
#' @param probs A numeric vector of probabilities between 0 and 1.
#' @return A scalar, vector, matrix, or nested list of quantiles.
#' @export
elementwise_quantile <- function(values, probs) {
  recursive_elementwise(values, quantile, probs = probs)
}

# -----------------------------------------------------------------------------
# Standard Error Computation Functions
# -----------------------------------------------------------------------------

#' Compute Standard Errors for Parameters in Transition Models
#'
#' Computes standard errors for various parameters in transition models
#' based on summaries extracted from a list of bootstrapped models.
#'
#' @param models A list of transition models from bootstrap.
#' @return A list of standard errors for each parameter.
#' @export
models_to_standard_errors <- function(models) {
  summaries <- summarize_transition_models(models)

  standard_errors <- list(
    treated_period = models[[1]]$treated_period,
    Y = models[[1]]$Y,
    prob_treated = elementwise_sd(lapply(models, `[[`, "prob_treated")),
    priors_treated = elementwise_sd(lapply(models, `[[`, "priors_treated")),
    pmfs_initial_treated = elementwise_sd(lapply(models, `[[`, "pmfs_initial_treated")),
    pmfs_initial_control = elementwise_sd(lapply(models, `[[`, "pmfs_initial_control")),
    pmfs_pre_treatment_treated = elementwise_sd(lapply(models, `[[`, "pmfs_pre_treatment_treated")),
    pmfs_pre_treatment_control = elementwise_sd(lapply(models, `[[`, "pmfs_pre_treatment_control")),
    ATT = elementwise_sd(lapply(summaries, `[[`, "ATT")),
    LTATTs = elementwise_sd(lapply(summaries, `[[`, "LTATTs")),
    Ps_treated = elementwise_sd(lapply(models, `[[`, "Ps_treated")),
    Ps_control = elementwise_sd(lapply(models, `[[`, "Ps_control")),
    Ps_treated_empirical = elementwise_sd(lapply(models, `[[`, "Ps_treated_empirical")),
    Ps_control_empirical = elementwise_sd(lapply(models, `[[`, "Ps_control_empirical"))
  )

  return(standard_errors)
}

#' Compute Quantiles for Parameters in Transition Models
#'
#' Computes quantiles for various parameters in transition models
#' based on summaries extracted from a list of bootstrapped models.
#'
#' @param models A list of transition models from bootstrap.
#' @param probs A numeric vector of probabilities between 0 and 1.
#' @return A list of quantiles for each parameter.
#' @export
models_to_quantiles <- function(models, probs) {
  summaries <- summarize_transition_models(models)

  quantiles <- list(
    treated_period = models[[1]]$treated_period,
    Y = models[[1]]$Y,
    prob_treated = elementwise_sd(lapply(models, `[[`, "prob_treated")),
    priors_treated = elementwise_sd(lapply(models, `[[`, "priors_treated")),
    pmfs_initial_treated = elementwise_sd(lapply(models, `[[`, "pmfs_initial_treated")),
    pmfs_initial_control = elementwise_sd(lapply(models, `[[`, "pmfs_initial_control")),
    pmfs_pre_treatment_treated = elementwise_sd(lapply(models, `[[`, "pmfs_pre_treatment_treated")),
    pmfs_pre_treatment_control = elementwise_sd(lapply(models, `[[`, "pmfs_pre_treatment_control")),
    ATT = elementwise_quantile(lapply(summaries, `[[`, "ATT"), probs = probs),
    LTATTs = elementwise_quantile(lapply(summaries, `[[`, "LTATTs"), probs = probs),
    Ps_treated = elementwise_quantile(lapply(models, `[[`, "Ps_treated"), probs = probs),
    Ps_control = elementwise_quantile(lapply(models, `[[`, "Ps_control"), probs = probs),
    Ps_treated_empirical = elementwise_quantile(lapply(models, `[[`, "Ps_treated_empirical"), probs = probs),
    Ps_control_empirical = elementwise_quantile(lapply(models, `[[`, "Ps_control_empirical"), probs = probs),
    Ps_diffs = elementwise_quantile(lapply(summaries, `[[`, "Ps_diffs"), probs = probs),
    Ps_empirical_diffs = elementwise_quantile(lapply(summaries, `[[`, "Ps_empirical_diffs"), probs = probs)
  )

  return(quantiles)
}
