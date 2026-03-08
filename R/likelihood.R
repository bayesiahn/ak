#' @title Estimate the log likelihood of a sequence of outcomes for a single pattern
#' @description
#'  This function estimates the log likelihood of a sequence of outcomes
#'  given a transition model of initial_outcome_dist and Ps_control
#'  for each group. The likelihood is given by
#'  \deqn{f_0^j(y_0) \prod_{t=1}^{T-1} f_t^j(y_t | y_{t-1})}
#' @param y_indices A vector of indices representing the sequence of outcomes.
#' @param initial_outcome_dist_j A K-length vector of JOINT initial outcome distributions (Y_1, D) for group j, fixing treatment status D
#' @param Ps_control_j A T-length list of transition matrices for group j.
#' @return The log likelihood of the sequence of outcomes.
#' @export
#' @examples
#' y_indices <- c(1, 2, 1)
#' initial_outcome_dist_j <- c(0.5, 0.5)
#' Ps_control_j <- list(matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2),
#'   matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2))
#' transition_to_log_likelihood_j(y_indices, initial_outcome_dist_j, Ps_control_j)
transition_to_log_likelihood_j <- function(y_indices, initial_outcome_dist_j, Ps_control_j) {
  transitions_to_log_likelihoods_j(matrix(y_indices, nrow = 1), initial_outcome_dist_j, Ps_control_j)
}

#' @title Estimate the log likelihoods of sequences of outcomes
#' @description
#' This function estimates the full log likelihood of a sequence of outcomes
#' given a transition model of initial_outcome_dist and Ps_control
#' for each group. The likelihood is given by
#' \deqn{\sum_{j=1}^J \pi_j f_0^j(y_0) \prod_{t=1}^{T-1} f_t^j(y_t | y_{t-1})}
#' @param y_indices_matrix N by (T+1) matrix of indices whose ith row represents ith unit's sequence of outcomes.
#' @param initial_outcome_dist_j A K-length vector of JOINT initial outcome distributions (Y_1, D) for group j, fixing treatment status D
#' @param Ps_control_j A T-length list of transition matrices for group j.
#' @return The log likelihood of the sequence of outcomes.
#' @export
#' @examples
#' y_indices_matrix <- matrix(c(1, 2, 1, 2, 2, 2), nrow = 2)
#' initial_outcome_dist_j <- c(0.5, 0.5)
#' Ps_control_j <- list(matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2),
#'  matrix(c(0.7, 0.5, 0.3, 0.5), nrow = 2))
#' transitions_to_log_likelihoods_j(y_indices_matrix, initial_outcome_dist_j, Ps_control_j)
transitions_to_log_likelihoods_j <- function(y_indices_matrix, initial_outcome_dist_j, Ps_control_j) {
  # Estimate likelihood
  log_likelihood <- log(initial_outcome_dist_j[y_indices_matrix[,1]])
  for (t in 2:ncol(y_indices_matrix)) {
    log_likelihood <- log_likelihood +
      log(Ps_control_j[[t-1]][y_indices_matrix[,(t-1):t,drop = FALSE]])
  }
  return(log_likelihood)
}

#' @title Estimate the log likelihood of a sequence of outcomes
#' @description
#'  This function estimates the full log likelihood of a sequence of outcomes
#'  given a transition model of initial_outcome_dist and Ps_control.
#' @param y_indices A vector of indices representing the sequence of outcomes.
#' @param initial_outcome_dists A J-length list of JOINT initial outcome distributions (Y_1, D) for each group, fixing treatment status D.
#' @param Ps_control A T-length list of transition matrices for group j.
#' @param priors A J-length vector of prior probabilities for each group.
#' @return The log likelihood of the sequence of outcomes.
#' @export
#' @examples
#' \dontrun{
#' y_indices <- c(1, 2, 1)
#' priors <- c(0.3, 0.7)
#' transition_model <- generate_transition_model(c(0,1), length(priors), length(y_indices)-1)
#' initial_outcome_dists <- transition_model$initial_outcome_dists
#' Ps_control <- transition_model$Ps_control
#' transition_to_log_likelihood(y_indices, initial_outcome_dists, Ps_control, priors)
#' }
transition_to_log_likelihood <- function(y_indices, initial_outcome_dists, Ps_control, priors) {
  return(transitions_to_log_likelihood(matrix(y_indices, nrow = 1),
                                  initial_outcome_dists, Ps_control, priors))
}

#' @title Estimate the log likelihood of sequences of outcomes
#' @description
#' This function estimates the full log likelihood of sequences of outcomes
#' given a transition model of initial_outcome_dist and Ps_control, optionally applying unit weights.
#' @param y_indices_matrix N by (T+1) matrix of indices whose ith row represents ith unit's sequence of outcomes.
#' @param initial_outcome_dists A J-length list of JOINT initial outcome distributions (Y_1, D) for each group, fixing treatment status D.
#' @param Ps_control A J-length list of transition matrices for each group.
#' @param priors A J-length vector of prior probabilities for each group.
#' @param unit_weights A numeric vector of weights for each unit, or \code{NULL}.
#'        If \code{NULL}, equal weights are assumed.
#' @return The log likelihood of the sequence of outcomes.
#' @export
#' @examples
#' \dontrun{
#' y_indices_matrix <- matrix(c(1, 2, 1, 2, 2, 2), nrow = 2)
#' priors <- c(0.3, 0.7)
#' transition_model <- generate_transition_model(c(0,1), length(priors), length(y_indices)-1)
#' initial_outcome_dists <- transition_model$initial_outcome_dists
#' Ps_control <- transition_model$Ps_control
#' transitions_to_log_likelihood(y_indices_matrix, initial_outcome_dists, Ps_control, priors)
#' }
transitions_to_log_likelihood <- function(y_indices_matrix, initial_outcome_dists, Ps_control, priors,
                                          unit_weights = NULL) {
  # Check for empty matrix
  if (length(y_indices_matrix) == 0) {
    return(0)
  }

  # Validate and handle unit_weights
  if (is.null(unit_weights)) {
    unit_weights <- rep(1, nrow(y_indices_matrix))
  }
  if (length(unit_weights) != nrow(y_indices_matrix)) {
    stop("Length of unit_weights must match the number of rows in y_indices_matrix.")
  }

  # Estimate likelihood for each latent group
  log_likelihood_js <- sapply(1:length(priors), function(j) {
    transitions_to_log_likelihoods_j(y_indices_matrix, initial_outcome_dists[[j]], Ps_control[[j]])
  })
  log_likelihood_js <- matrix(log_likelihood_js, nrow = nrow(y_indices_matrix))

  # Estimate full likelihood, applying unit weights
  log_likelihoods <- apply(log_likelihood_js, 1, function(log_likelihood_j) {
    log(sum(priors * exp(log_likelihood_j)))
  })

  # Apply weights to log likelihoods
  weighted_log_likelihood <- sum(log_likelihoods * unit_weights)
  return(weighted_log_likelihood)
}

#' @title Estimate the posterior probability of  sequence of outcomes
#' @description
#' This function estimates the posterior probability of a sequence of outcomes
#' given a transition model of initial_outcome_dist and Ps_control
#' for each group.
#' @param y_indices_matrix N by (T+1) matrix of indices whose ith row represents ith unit's sequence of outcomes.
#' @param initial_outcome_dists A J-length list of JOINT initial outcome distributions (Y_1, D) for each group, fixing treatment status D.
#' @param Ps A J-length list of transition matrices for each group.
#' @param priors A J-length vector of prior probabilities for each group.
#' @return A N by J matrix of posterior probabilities for each individual belonging to each group.
#' @export
#' @examples
#' \dontrun{
#' priors <- c(0.3, 0.7)
#' transition_model <- generate_transition_model(c(0,1), length(priors), length(y_indices)-1)
#' initial_outcome_dists <- transition_model$initial_outcome_dists
#' Ps_control <- transition_model$Ps_control
#' transitions_to_posteriors(y_indices_matrix, initial_outcome_dists, Ps_control, priors)
#' }
transitions_to_posteriors <- function(y_indices_matrix, initial_outcome_dists, Ps, priors) {
  J <- length(priors)

  # If no observations, return empty matrix
  if (length(y_indices_matrix) == 0) {
    return(matrix(0, nrow = 0, ncol = J))
  }

  # If only one component, return matrix of ones
  if (J == 1) {
    return(matrix(1, nrow = nrow(y_indices_matrix), ncol = 1))
  }

  # Compute log lik for each individual i and transition pattern j; N by J matrix
  log_lik_js <- sapply(1:J, function(j)
    transitions_to_log_likelihoods_j(y_indices_matrix, initial_outcome_dists[[j]], Ps[[j]]))
  log_lik_js <- matrix(log_lik_js, ncol = J)

  # For each individual, find the maximum log-likelihood across all j components
  log_lik_max <- apply(log_lik_js, 1, max)

  # Compute numerator for each ith individual
  posterior_numerators <- sapply(1:J, function(j)
    priors[j] * exp(log_lik_js[,j] - log_lik_max))
  posterior_numerators <- matrix(posterior_numerators, ncol = J)

  # Compute posterior probability of each i belonging to each j compaonent
  posterior_numerators / apply(posterior_numerators, 1, sum)
}
