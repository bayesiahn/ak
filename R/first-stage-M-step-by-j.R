
#' Compute Initial Distribution for a Specific Latent Group
#'
#' Given posterior probabilities and initial outcome indicators, this function computes
#' the initial distribution for a specific latent group. It calculates the weighted sum
#' of the initial outcome indicators and normalizes them to form a probability distribution.
#'
#' @param posteriors_j A vector of posterior probabilities for a specific latent group.
#' @param initial_outcome_indicators A matrix where each row is a one-hot encoded
#'        vector representing the initial state of an individual.
#'
#' @return A normalized vector representing the initial distribution for the specified
#'         latent group.
#'
#' @examples
#' \dontrun{
#' posteriors_j <- runif(2)
#' initial_outcome_indicators <- matrix(c(1, 0, 0, 1), nrow = 2)
#' initial_outcome_dists_j <- M_step_initial_outcome_dists_j(posteriors_j, initial_outcome_indicators)
#' }
M_step_initial_outcome_dists_j <- function(posteriors_j, initial_outcome_indicators) {
  # Check for empty posteriors
  if (length(posteriors_j) == 0) {
    return(numeric(0))
  }

  # Compute weighted sum for each outcome using efficient O(n) operation
  # crossprod(x, y) = t(x) %*% y, so t(posteriors) %*% indicators gives 1×K weighted sum
  # as.matrix() ensures sparse matrices (from Matrix::bdiag) are converted to dense
  initial_outcome_dists_j <- as.vector(crossprod(posteriors_j, as.matrix(initial_outcome_indicators)))

  # Normalize and return
  return(initial_outcome_dists_j / sum(initial_outcome_dists_j))
}

#' Compute Transition Matrices for Each Period for a Specific Latent Group
#'
#' Given posterior probabilities and transition indicators, this function computes
#' the transition matrices for each period for a specific latent group.
#'
#' @param posteriors_j A vector of posterior probabilities for a specific latent group.
#' @param transition_indicators A list of transition indicators for each time period.
#'
#' @return A list of transition matrices for each period for the specified latent group.
#'
#' @examples
#' \dontrun{
#' posteriors_j <- runif(2)
#' transition_indicators <- list(matrix(1:4, nrow = 2), matrix(5:8, nrow = 2))
#' Ps_j <- M_step_Ps_j(posteriors_j, transition_indicators)
#' }
M_step_Ps_j <- function(posteriors_j, transition_indicators) {
  # Compute T_max from transition_indicators
  T_max <- length(transition_indicators)

  # Compute transition matrices for each period
  return(lapply(1:T_max, function(t) {
    M_step_P_jt(posteriors_j, transition_indicators[[t]])
  }))
}

#' Compute Transition Matrix for a Specific Period and Latent Group
#'
#' This function computes the transition matrix for a specific time period
#' and latent group, using the provided posterior probabilities and transition indicators.
#'
#' @param posteriors_j A vector of posterior probabilities for a specific latent group.
#' @param transition_indicators_t Transition indicators for a specific time period.
#'
#' @return A transition matrix for the specified time period and latent group.
#'
#' @examples
#' \dontrun{
#' posteriors_j <- matrix(runif(4), ncol = 2)
#' transition_indicators_t <- matrix(1:4, nrow = 2)
#' P_jt <- M_step_P_jt(posteriors_j, transition_indicators_t)
#' }
M_step_P_jt <- function(posteriors_j, transition_indicators_t) {
  # Step 1: Multiply each matrix in transition_indicator by posterior weights
  P_jt <- mapply(function(transition_indicator, posterior_weight)
    transition_indicator * posterior_weight,
    transition_indicators_t, posteriors_j, SIMPLIFY = FALSE)

  # Step 2: Sum the matrices element-wise
  P_jt <- Reduce(`+`, P_jt)

  # Step 3: Normalize each row of the matrix
  # If row sums to zero (no transitions observed from that state),
  # use uniform distribution to maintain valid stochastic matrix
  n_outcomes <- ncol(P_jt)
  P_jt <- t(apply(P_jt, 1, function(row) {
    row_sum <- sum(row)
    if (row_sum == 0) {
      # No transitions observed: use uniform distribution
      return(rep(1 / n_outcomes, n_outcomes))
    } else {
      return(row / row_sum)
    }
  }))

  return(P_jt)
}
