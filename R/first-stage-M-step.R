
#' Perform the M-Step of an EM Algorithm
#'
#' Computes priors, initial outcome distributions, transition matrices, and
#' the complete-data log-likelihood for the current parameter estimates
#' given the posterior probabilities.
#'
#' @param posteriors A matrix of posterior probabilities (n × J), where each row
#'        corresponds to a unit and each column to a latent group.
#' @param data_for_est A list returned by \code{split_data_for_estimation()}, containing
#'        state sequences and indicators split by treatment status.
#' @param unit_weights Optional numeric vector of weights for each unit. If \code{NULL},
#'        equal weights are used.
#'
#' @return A list containing:
#' \describe{
#'   \item{priors}{Estimated prior probabilities over latent groups.}
#'   \item{priors_treated}{Estimated prior probabilities over latent groups for treated units.}
#'   \item{pmfs_initial_control}{Initial outcome distributions for control units.}
#'   \item{pmfs_initial_treated}{Initial outcome distributions for treated units.}
#'   \item{Ps_control}{Estimated transition matrices for control units.}
#'   \item{Ps_treated}{Estimated transition matrices for treated units.}
#'   \item{log_likelihood}{Complete-data log-likelihood.}
#' }
#'
#' @export
M_step <- function(posteriors, data_for_est, unit_weights = NULL) {
  unit_weights <- check_unit_weights_validity(posteriors, unit_weights)

  # Compute number of latent groups and fetch indicators
  J <- ncol(posteriors)

  is_control <- data_for_est$is_control
  is_treated <- data_for_est$is_treated
  weights_control <- unit_weights[is_control]
  weights_treated <- unit_weights[is_treated]

  # Split posteriors and weights
  posteriors_control <- matrix(posteriors[is_control, ], ncol = J)
  posteriors_treated <- matrix(posteriors[is_treated, ], ncol = J)

  # Estimate priors
  priors <- M_step_priors(posteriors, unit_weights)
  priors_treated <- M_step_priors(posteriors_treated, weights_treated)

  # Estimate initial outcome distributions
  pmfs_initial_control <- M_step_initial_outcome_dists(
    posteriors_control,
    data_for_est$initial_outcome_indicators_control,
    weights_control
  )
  pmfs_initial_treated <- M_step_initial_outcome_dists(
    posteriors_treated,
    data_for_est$initial_outcome_indicators_treated,
    weights_treated
  )

  # Estimate transition matrices (shared in pretreatment periods)
  Ps_pre <- M_step_Ps(posteriors, data_for_est$transition_indicators_pretreatment, unit_weights)

  # Estimate transition matrices for control and treated units in post-treatment periods
  Ps_control_post <- M_step_Ps(posteriors_control, data_for_est$transition_indicators_control_post, weights_control)
  Ps_treated_post <- M_step_Ps(posteriors_treated, data_for_est$transition_indicators_treated_post, weights_treated)

  # Concatenate transition matrices; note that Ps_control and Ps_treated are the same in pretreatment periods
  Ps_control <- concatenate_Ps(Ps_pre, Ps_control_post)
  Ps_treated <- concatenate_Ps(Ps_pre, Ps_treated_post)

  # Initial P(Y_1, D | Z) distribution
  p_y1d <- M_step_initial_outcome_dists(
    rbind(posteriors_control, posteriors_treated),
    data_for_est$initial_outcome_indicators_joint,
    c(weights_control, weights_treated)
  )
  # Make it into a matrix form
  p_y1d <- lapply(p_y1d, function(p_y1d_j) matrix(p_y1d_j, ncol = 2))
  p_y1d_splitted <- split_p_y1d(p_y1d)

  # Log-likelihood computation
  ll_control <- transitions_to_log_likelihood(
    data_for_est$y_indices_control,
    p_y1d_splitted$p_y1d_control,
    Ps_control,
    priors,
    weights_control
  )
  ll_treated <- transitions_to_log_likelihood(
    data_for_est$y_indices_treated,
    p_y1d_splitted$p_y1d_treated,
    Ps_treated,
    priors,
    weights_treated
  )

  log_likelihood <- sum(ll_control) + sum(ll_treated)

  list(
    priors = priors,
    priors_treated = priors_treated,
    pmfs_initial_control = pmfs_initial_control,
    pmfs_initial_treated = pmfs_initial_treated,
    Ps_control = Ps_control,
    Ps_treated = Ps_treated,
    p_y1d = p_y1d,
    log_likelihood = log_likelihood
  )
}

#' Compute Prior Probabilities for Each Latent Group
#'
#' This function calculates the prior probabilities for each latent group based on
#' the provided posterior probabilities. It normalizes the sum of posteriors across
#' all groups to obtain the priors, optionally applying unit weights.
#'
#' @param posteriors A matrix of posterior probabilities where each column
#'        corresponds to a latent group.
#' @param unit_weights A numeric vector of weights for each unit, with the same length
#'        as the number of rows in \code{posteriors}. If \code{NULL}, equal weights are assumed.
#'
#' @return A numeric vector representing the prior probabilities for each latent group.
#' @examples
#' posteriors <- matrix(c(0.2, 0.3, 0.4, 0.1), ncol = 2)
#' priors <- M_step_priors(posteriors)
#' print(priors)  # Should return normalized priors for each group
#'
#' unit_weights <- generate_weights(nrow(posteriors))
#' priors_weighted <- M_step_priors(posteriors, unit_weights = unit_weights)
#' print(priors_weighted)  # Should return priors accounting for unit weights
#' @export
M_step_priors <- function(posteriors, unit_weights = NULL) {
  # Validate and handle unit_weights
  unit_weights <- check_unit_weights_validity(posteriors, unit_weights)

  # Compute weighted sum for each column
  priors <- colSums(posteriors * unit_weights)

  # Normalize and return
  return(priors / sum(priors))
}

#' Compute Initial Distributions for Each Latent Group
#'
#' This function calculates the initial distributions for each latent group
#' based on posteriors and initial outcome indicators. The function applies
#' `M_step_initial_outcome_dists_j` across all latent groups, optionally using
#' unit weights for weighted calculations.
#'
#' @param posteriors A matrix of posterior probabilities for each individual
#'        belonging to each group.
#' @param initial_outcome_indicators A list of one-hot encoded vectors representing
#'        the initial state of each individual.
#' @param unit_weights A numeric vector of weights for each unit, or \code{NULL}.
#'        If \code{NULL}, equal weights are assumed.
#'
#' @return A list of initial distributions, one for each latent group.
#' @examples
#' posteriors <- matrix(runif(6), nrow = 2)
#' initial_outcome_indicators <- rbind(c(1, 0), c(0, 1))
#' initial_dists <- M_step_initial_outcome_dists(posteriors, initial_outcome_indicators)
#' print(initial_dists)
#' @export
M_step_initial_outcome_dists <- function(posteriors, initial_outcome_indicators,
                                         unit_weights = NULL) {
  # Validate and handle unit weights
  unit_weights <- check_unit_weights_validity(posteriors, unit_weights)

  # Check for empty posteriors
  if (length(initial_outcome_indicators) < 2) {
    if (length(initial_outcome_indicators) == 0) {
      return(NULL)
    }
    posteriors <- matrix(posteriors, nrow = 1)
  }

  # Compute J from initial_outcome_indicators
  posteriors <- as.matrix(posteriors)
  J <- ncol(posteriors)

  # Compute initial distributions for each transition pattern, incorporating unit weights
  return(lapply(1:J, function(j) {
    M_step_initial_outcome_dists_j(posteriors[, j] * unit_weights, initial_outcome_indicators)
  }))
}

#' Compute Transition Matrices for Each Transition Pattern
#'
#' This function computes transition matrices for each latent group based on
#' posteriors, transition indicators, and optional unit weights.
#'
#' @param posteriors A matrix of posterior probabilities for each individual belonging to each group.
#' @param transition_indicators A list of transition indicators for each time period.
#' @param unit_weights A numeric vector of weights for each unit, or \code{NULL}.
#'        If \code{NULL}, equal weights are assumed.
#'
#' @return A list of transition matrices for each latent group.
#' @examples
#' posteriors <- matrix(runif(6), nrow = 2)
#' transition_indicators <- list(matrix(1:4, nrow = 2), matrix(5:8, nrow = 2))
#' unit_weights <- c(1, 2)
#' Ps <- M_step_Ps(posteriors, transition_indicators, unit_weights)
#' print(Ps)
#' @export
M_step_Ps <- function(posteriors, transition_indicators, unit_weights = NULL) {
  # Validate and handle unit weights
  unit_weights <- check_unit_weights_validity(posteriors, unit_weights)

  # Compute J from posteriors
  J <- ncol(posteriors)

  # Compute transition matrices for each latent group, applying unit weights
  return(lapply(1:J, function(j) {
    M_step_Ps_j(posteriors[, j] * unit_weights, transition_indicators)
  }))
}
