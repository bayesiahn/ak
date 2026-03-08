#' Convert Score Matrix to Probability Mass Function (PMF)
#'
#' This function converts an \eqn{N \times (K-\mathrm{minimal})} matrix of score values into
#' a probability mass function (PMF) by normalizing across observations,
#' adjusting it based on an initial PMF.
#'
#' @param pmf_j_score A numeric \eqn{N \times (K-\mathrm{minimal})} matrix of score values, where
#'   each row corresponds to an observation and each column corresponds to
#'   a possible outcome category.
#' @param pmf_j A numeric vector of length \eqn{K}, representing the initial PMF.
#' @param weights An optional numeric vector of length \eqn{N}, representing
#'   observation weights. If provided, scores are weighted accordingly.
#' @param minimal A logical value. If TRUE, expects \eqn{N \times (K-1)} matrix
#' for scores.
#'
#' @return A numeric vector of length \eqn{K}, representing the updated PMF.
#'
#' @examples
#' set.seed(123)
#' N <- 100
#' K <- 3
#' pmf_j_score <- matrix(runif(N * K), nrow = N, ncol = K)
#' pmf_j <- rep(1/K, K)  # Uniform initial PMF
#' weights <- runif(N)
#'
#' pmf <- score_to_pmf(pmf_j_score, pmf_j, weights)
#' print(sum(pmf))  # Should be a normalized PMF
score_to_pmf <- function(pmf_j_score, pmf_j, weights = NULL, minimal = TRUE) {
  # Compute the number of units
  N <- nrow(pmf_j_score)

  # If weights are not NULL, multiply each row of scores by the weights
  if (!is.null(weights)) {
    if (length(weights) != N) {
      stop("Length of weights must be equal to the number of rows in pmf_j_score")
    }
    pmf_j_score <- pmf_j_score * weights
  }

  # Return the updated PMF
  # If minimal option is TRUE, the last element is set to be the 1-sum(rest)
  score_sums <- colSums(pmf_j_score)
  if (minimal) {
    score_sums <- c(score_sums, N - sum(score_sums))
  }
  pmf <- pmf_j + (score_sums / N)

  # Normalize
  pmf <- pmax(pmf, 0)
  pmf / sum(pmf)
}

#' Convert List of Score Matrices to Probability Mass Functions (PMFs)
#'
#' This function converts a list of score matrices into a list of updated
#' probability mass functions (PMFs) by normalizing across observations.
#'
#' @param pmfs_score A list of \eqn{J} elements, where each element is an \eqn{N \times (K-\mathrm{minimal})}
#'   matrix of score values, with each row corresponding to an observation and each column
#'   corresponding to a possible outcome category.
#' @param pmfs A list of \eqn{J} numeric vectors, each of length \eqn{K}, representing the
#'   initial PMFs.
#' @param weights An optional numeric vector of length \eqn{N}, representing observation weights.
#'   If provided, scores are weighted accordingly.
#' @param minimal A logical value. If TRUE, expects \eqn{N \times (K-1)} matrix
#' for scores.
#'
#' @return A list of \eqn{J} numeric vectors, each of length \eqn{K}, representing the updated PMFs.
#'
#' @examples
#' set.seed(123)
#' N <- 100
#' K <- 3
#' J <- 2
#' pmfs_score <- lapply(1:J, function(x) matrix(runif(N * K), nrow = N, ncol = K))
#' pmfs <- lapply(1:J, function(x) rep(1/K, K))  # Uniform initial PMFs
#' weights <- runif(N)
#'
#' updated_pmfs <- score_to_pmfs(pmfs_score, pmfs, weights)
#' print(length(updated_pmfs))  # Should be J
#' print(length(updated_pmfs[[1]]))  # Should be K
score_to_pmfs <- function(pmfs_score, pmfs, weights = NULL, minimal = TRUE) {
  J <- length(pmfs_score)
  lapply(1:J, function(j) score_to_pmf(as.matrix(pmfs_score[[j]]), pmfs[[j]], weights, minimal))
}

#' Convert Score List to Transition Probability Matrix
#'
#' This function converts a list of \eqn{N \times K_{t-1}} score matrices into
#' a \eqn{K_t \times K_{t-1}} transition probability matrix, adjusting it based on
#' an initial transition probability matrix.
#'
#' @param P_jt_score A list of \eqn{K_t} elements, where each element is an
#'   \eqn{N \times (K_{t-1}-\mathrm{minimal})} matrix of score values.
#' @param P_jt A \eqn{K_t \times K_{t-1}} matrix representing the initial transition probabilities.
#' @param weights An optional numeric vector of length \eqn{N}, representing
#'   observation weights. If provided, scores are weighted accordingly.
#' @param minimal A logical value. If TRUE, expects \eqn{N \times (K_{t-1}-1)} matrix
#' for scores.
#'
#' @return A \eqn{K_t \times K_{t-1}} transition probability matrix.
#'
#' @examples
#' set.seed(123)
#' N <- 100
#' K_t <- 3
#' K_t_minus_1 <- 3
#' P_jt_score <- lapply(1:K_t, function(x) matrix(runif(N * K_t_minus_1), nrow = N, ncol = K_t_minus_1))
#' P_jt <- matrix(1/K_t_minus_1, K_t, K_t_minus_1)  # Uniform initial transitions
#' weights <- runif(N)
#'
#' updated_P_jt <- score_to_P_jt(P_jt_score, P_jt, weights)
#' print(dim(updated_P_jt))  # Should be (K_t, K_{t-1})
score_to_P_jt <- function(P_jt_score, P_jt, weights = NULL,  minimal = TRUE) {
  # Compute the length of past outcome support
  K_t <- length(P_jt_score)

  # Convert each row of P_jt_score to a transition probability distribution
  P_jt <- sapply(1:K_t, function(k_t) {
    score_to_pmf(P_jt_score[[k_t]], P_jt[k_t,], weights, minimal)
  })

  # Return as a matrix; transpose needed because sapply returns columns
  t(P_jt)
}

#' Convert Score List to Time-Series Transition Matrices
#'
#' This function converts a list of time-step transition score matrices into
#' a time-series transition probability model, adjusting it based on an initial
#' set of transition probability matrices.
#'
#' @param P_j_score A list of \eqn{T-1} elements, where each element is a list of
#'   \eqn{K_t} matrices representing score values for transitions at time step \eqn{t}.
#' @param P_j A list of \eqn{T-1} elements, where each element is a \eqn{K_t \times K_{t-1}}
#'   transition probability matrix representing the initial transition probabilities at time step \eqn{t}.
#' @param weights An optional numeric vector of length \eqn{N}, representing
#'   observation weights. If provided, scores are weighted accordingly.
#' @param minimal A logical value. If TRUE, expects \eqn{N \times (K-1)} matrix
#' for scores.
#'
#' @return A list of \eqn{T-1} transition probability matrices, where each matrix represents
#'   the updated transition probabilities for a specific time step.
#'
#' @examples
#' set.seed(123)
#' T_max <- 3
#' K_t <- 3
#' K_t_minus_1 <- 3
#' P_j_score <- lapply(1:(T_max-1), function(x)
#'   lapply(1:K_t, function(y) matrix(runif(100 * K_t_minus_1), nrow = 100, ncol = K_t_minus_1))
#' )
#' P_j <- lapply(1:(T_max-1), function(x) matrix(1/K_t_minus_1, K_t, K_t_minus_1))  # Uniform initial transitions
#' weights <- runif(100)
#'
#' updated_P_j <- score_to_P_j(P_j_score, P_j, weights)
#' print(length(updated_P_j))  # Should be T-1
score_to_P_j <- function(P_j_score, P_j, weights = NULL, minimal = TRUE) {
  T_max <- length(P_j_score) + 1
  lapply(1:(T_max-1), function(t) score_to_P_jt(P_j_score[[t]], P_j[[t]], weights, minimal))
}

#' Convert Score List to Transition Models for All Components
#'
#' This function converts a list of transition score values across multiple mixture
#' components into a time-series transition probability model.
#'
#' @param Ps_score A list of \eqn{J} elements, where each element is a list of
#'   \eqn{T-1} transition probability score matrices.
#' @param Ps A list of \eqn{J} elements, where each element is a list of \eqn{T-1} transition
#'   probability matrices representing the initial transitions.
#' @param minimal A logical value. If TRUE, expects \eqn{N \times (K-1)} matrix
#' for scores.
#' @param weights An optional numeric vector of length \eqn{N}, representing
#'   observation weights. If provided, scores are weighted accordingly.
#'
#' @return A list of \eqn{J} transition probability models.
#'
#' @examples
#' set.seed(123)
#' T_max <- 3
#' J <- 2
#' K_t <- 3
#' K_t_minus_1 <- 3
#' Ps_score <- lapply(1:J, function(j)
#'   lapply(1:(T_max-1), function(x)
#'     lapply(1:K_t, function(y) matrix(runif(100 * K_t_minus_1), nrow = 100, ncol = K_t_minus_1))
#'   )
#' )
#' Ps <- lapply(1:J, function(j)
#'   lapply(1:(T_max-1), function(x) matrix(1/K_t_minus_1, K_t, K_t_minus_1))
#' )
#' weights <- runif(100)
#'
#' updated_Ps <- score_to_Ps(Ps_score, Ps, weights)
#' print(length(updated_Ps))  # Should be J
score_to_Ps <- function(Ps_score, Ps, weights = NULL, minimal = TRUE) {
  J <- length(Ps_score)
  lapply(1:J, function(j) score_to_P_j(Ps_score[[j]], Ps[[j]], weights, minimal))
}
