#' Compute Score Function for Initial Outcome Distribution
#'
#' This function calculates the score function for the initial outcome
#' distribution of a mixture component \eqn{j}.
#'
#' @param y_indices Integer vector of length \eqn{N}, where each entry
#'   represents the observed category index (values in \eqn{\{1, \ldots, K\}}).
#' @param pmf_j Numeric vector of length \eqn{K}
#'   representing the initial probability distribution of outcomes.
#' @param posteriors_j Numeric vector of length \eqn{N}, representing
#'   posterior probabilities for component \eqn{j} for each observation.
#' @param minimal A logical value. If TRUE, returns \eqn{N \times (K-1)} matrix.
#'
#' @return A numeric \eqn{N \times (K-\mathrm{minimal})} matrix, where each row corresponds
#'   to an observation and each column represents the score contribution
#'   for each outcome category.
#'
#' @examples
#' y_indices <- c(1, 2, 1, 3)
#' pmf_j <- c(0.2, 0.5, 0.3)
#' posteriors_j <- c(0.7, 0.8, 0.6, 0.9)
#' score_pmf_j(y_indices, pmf_j, posteriors_j)
score_pmf_j <- function(y_indices, pmf_j, posteriors_j, minimal = TRUE) {
  # Number of possible outcome categories (K)
  K <- length(pmf_j)

  # Number of observations (N)
  N <- length(y_indices)

  # Check if K = 1 (single category)
  if (K <= 1) {
    # Return a matrix of zeros (no variation)
    return(matrix(0, nrow = N, ncol = K-minimal))
  }

  # Initialize a matrix to store score values (N x K)
  score_j <- matrix(0, nrow = N, ncol = K)

  # Assign probabilities from pmf_j to score_j at the observed indices
  p_minimum <- min(min(c(1, pmf_j[pmf_j != 0])), 1/N)  # Minimum non-zero probability
  pmf_j_enforced <- pmax(pmf_j, p_minimum)  # Enforce non-zero probabilities
  score_j[cbind(1:N, y_indices)] <- 1/pmf_j_enforced[y_indices]

  # Adjust the score matrix: subtract the last column (K) from each column
  # Use vectorized operations instead of apply() for efficiency
  score_j <- (score_j - score_j[, K]) * posteriors_j

  # Ensure last column (K) accounts for remaining probability mass
  score_j[, K] <- 1 - rowSums(matrix(score_j[, 1:(K-1)], ncol = (K-1)))

  # Choose appropriate columns based on minimal flag
  score_j <- matrix(score_j[, 1:(K-minimal)], ncol = (K-minimal))

  # Set colnames
  colnames(score_j) <- paste0("y", c(1:(K-minimal)))

  # Return the computed score function matrix
  return(score_j)
}

#' Compute Score Function for Mixture Component PMFs
#'
#' This function calculates the score function for the probability mass functions (PMFs)
#' of multiple mixture components, given observed outcome indices.
#'
#' @param y_indices Integer vector of length \eqn{N}, where each entry
#'   represents the observed category index (values in \eqn{\{1, \ldots, K\}}).
#' @param pmfs A list of \eqn{J} numeric vectors, where each vector has length \eqn{K},
#'   representing the initial probability mass function (PMF) for each mixture component.
#' @param posteriors A numeric \eqn{N \times J} matrix, where each entry represents
#'   the posterior probability of observation \eqn{i} belonging to component \eqn{j}.
#' @param t_index An optional time index to append to the column names.
#' @param minimal A logical value. If TRUE, returns \eqn{N \times (K-1)} matrices.
#'
#' @return A list of \eqn{J} matrices, where each matrix has dimensions \eqn{N \times (K-\mathrm{minimal})},
#'   representing the score contributions for each mixture component.
#'
#' @examples
#' y_indices <- c(1, 2, 1, 3)
#' pmfs <- list(c(0.2, 0.5, 0.3), c(0.3, 0.4, 0.3))
#' posteriors <- matrix(c(0.7, 0.8, 0.6, 0.9, 0.3, 0.2, 0.4, 0.1), nrow = 4, byrow = FALSE)
#'
#' scores <- score_pmfs(y_indices, pmfs, posteriors)
#' print(dim(scores[[1]]))  # Should be (4, 3)
#'
#' @export
score_pmfs <- function(y_indices, pmfs, posteriors, t_index = NULL, minimal = TRUE) {
  # Number of mixture components (J)
  J <- length(pmfs)

  # Initialize a list to store score matrices
  scores <- list()

  # Iterate over all mixture components
  for (j in 1:J) {
    # Compute the score matrix for the j-th component
    score_j <- score_pmf_j(y_indices, pmfs[[j]], posteriors[, j], minimal)

    # Add a time index to the column names
    if (ncol(score_j) > 0) {
      colnames(score_j) <- paste0("j", j, "_", colnames(score_j))
      if (!is.null(t_index)) {
        colnames(score_j) <- paste0("t", t_index, "_", colnames(score_j))
      }
    }

    # Store the score matrix for the j-th component
    scores[[j]] <- score_j
  }

  # Return the list of score matrices for all mixture components
  return(scores)
}

#' Compute Score Function for Transition Probabilities
#'
#' This function calculates the score function for a given transition
#' probability matrix \eqn{P_{jt}}, which represents the probabilities
#' of transitioning from \eqn{y_{t-1}} to \eqn{y_t}.
#'
#' @param y_indices_matrix An \eqn{N \times 2} matrix, where each row contains
#'   a past outcome \eqn{y_{t-1}} (first column) and a future outcome \eqn{y_t}
#'   (second column).
#' @param P_jt A \eqn{K \times K} matrix of transition probabilities,
#'   where \eqn{P_{jt}[k, m]} represents the probability of transitioning
#'   from \eqn{y_{t-1} = k} to \eqn{y_t = m}.
#' @param posteriors_j A numeric vector of length \eqn{N}, containing
#'   posterior probabilities for each observation.
#' @param j_index An optional group index to append to the column names.
#' @param t_index An optional time index to append to the column names.
#' @param minimal A logical value. If TRUE, returns \eqn{N \times (K-1)} matrices.
#'
#' @return A list of \eqn{K} matrices, where each element corresponds to
#'   a past outcome category \eqn{y_{t-1} = k}, and each matrix is of
#'   dimension \eqn{N \times (K-\mathrm{minimal})}, representing the score contributions for
#'   each possible transition \eqn{y_{t-1} \to y_t}.
#'
#' @examples
#' y_indices_matrix <- matrix(c(1, 2, 2, 3, 1, 1, 3, 2), ncol = 2, byrow = TRUE)
#' P_jt <- matrix(c(0.5, 0.5, 0,  0.3, 0.4, 0.3,  0.2, 0.3, 0.5), ncol = 3)
#' posteriors_j <- c(0.6, 0.8, 0.7, 0.9)
#'
#' scores <- score_P_jt(y_indices_matrix, P_jt, posteriors_j)
#' print(dim(scores[[1]]))  # Should be (N, K)
#'
#' @export
score_P_jt <- function(y_indices_matrix, P_jt, posteriors_j, j_index = NULL, t_index = NULL, minimal = TRUE) {
  # Check if y_indices_matrix has two columns
  if (ncol(y_indices_matrix) != 2) {
    stop("y_indices_matrix must have two columns to compute scores for P_jt")
  }

  # Number of possible outcome categories in the past outcomes
  K <- ncol(P_jt)

  # Number of observations (N); each row corresponds to an observation
  N <- nrow(y_indices_matrix)

  # Initialize a list to store score matrices
  scores_j <- list()

  # Iterate over all y_{t-1}
  for (k in 1:K) {
    # Find indices of y_{t-1} = k
    unit_indices_k <- which(y_indices_matrix[, 1] == k)

    # Initialize a matrix to store score values (N x (K-minimal))
    scores_jk <- matrix(0, nrow = N, ncol = (K-minimal))
    if ((length(P_jt[k,])-minimal) > 0) {
      colnames(scores_jk) <- paste0(1:(length(P_jt[k,])-minimal), "_", k)
    }

    # Check if there are any observations with y_{t-1} = k
    if (length(unit_indices_k) > 0) {
      # Find indices of y_t = 1, ..., K
      y_indices_k_future <- matrix(y_indices_matrix[unit_indices_k, 2])

      # Find probabilities of y_t = 1, ..., K given y_{t-1} = k
      P_jt_k <- P_jt[k, ]

      # Compute the score for each y_t
      scores_jk[unit_indices_k,] <- score_pmf_j(y_indices_k_future, P_jt_k, posteriors_j[unit_indices_k], minimal)
    }

    if (!is.null(t_index)) {
      colnames(scores_jk) <- paste0("t", t_index, "_", colnames(scores_jk))
    }
    if (!is.null(j_index)) {
      colnames(scores_jk) <- paste0("j", j_index, "_", colnames(scores_jk))
    }

    # Store the scores for y_{t-1} = k
    scores_j[[k]] <- scores_jk
  }

  # Return the scores for all y_{t-1} = 1, ..., K
  scores_j
}

#' Compute Score Function for All Transition Probabilities Across Time
#'
#' This function calculates the score function for all transition probability matrices
#' across a sequence of time steps, using individual time-step transition probabilities.
#'
#' @param y_indices_matrix An \eqn{N \times T} matrix, where each row represents
#'   an observation's sequence of outcomes across \eqn{T} time steps.
#' @param P_j A list of \eqn{T - 1} transition probability matrices, where each \eqn{P_j[[t]]}
#'   is a \eqn{K \times K} matrix containing transition probabilities at time step \eqn{t}.
#' @param posteriors_j A numeric vector of length \eqn{N}, containing posterior probabilities
#'   for each observation.
#' @param j_index An optional group index to append to the column names.
#' @param minimal A logical value. If TRUE, returns \eqn{N \times (K-1)} matrices.
#'
#' @return A list of \eqn{T - 1} elements, where each element is a list of \eqn{K} matrices
#'   (one for each past outcome category \eqn{y_{t-1}}), with dimensions \eqn{N \times (K-\mathrm{minimal})}.
#'   These matrices represent the score contributions for each transition \eqn{y_{t-1} \to y_t}.
#'
#' @examples
#' outcomes <- 1:2
#' T_max <- 3
#' model <- generate_transition_model(outcomes = outcomes, J = 1, T_max = T_max)
#' sample <- generate_sample(model, N = 100)
#' P_j <- model$Ps_control[[1]]  # List of transition matrices
#' y_indices_matrix <- sample$y
#' posteriors_j <- rep(1, nrow(y_indices_matrix))
#'
#' scores <- score_P_j(y_indices_matrix, P_j, posteriors_j)
#' print(length(scores))  # Should be T-1
#' print(dim(scores[[1]][[1]]))  # Should be (N, K)
#'
#' @export
score_P_j <- function(y_indices_matrix, P_j, posteriors_j, j_index = NULL, minimal = TRUE) {
  # Check if length(P_j) = ncol(y_indices_matrix) - 1
  if (length(P_j) != (ncol(y_indices_matrix) - 1)) {
    stop("P_j must have length equal to the number of columns in y_indices_matrix minus 1")
  }

  # Initialize a list to store score matrices for each time step
  scores_j <- list()

  # Iterate over all time steps
  for (t in 1:length(P_j)) {
    scores_j[[t]] <- score_P_jt(matrix(y_indices_matrix[, t:(t+1)], ncol = 2),
                                P_j[[t]], posteriors_j, j_index, t, minimal)
  }

  # Return the list of list of score matrices for all time steps
  scores_j
}

#' Compute Score Function for Transition Probabilities Across All Components
#'
#' This function calculates the score function for a set of transition probability
#' matrices across multiple mixture components, given sequences of observed outcomes.
#'
#' @param y_indices_matrix An \eqn{N \times T} matrix, where each row represents
#'   an observation's sequence of outcomes across \eqn{T} time steps.
#' @param Ps A list of \eqn{J} elements, where each \eqn{P_j} is a list of \eqn{T - 1}
#'   transition probability matrices (\eqn{K \times K}) corresponding to mixture component \eqn{j}.
#' @param posteriors A numeric \eqn{N \times J} matrix, where each entry represents
#'   the posterior probability of observation \eqn{i} belonging to component \eqn{j}.
#' @param minimal A logical value. If TRUE, returns \eqn{N \times (K-1)} matrices.
#'
#' @return A list of \eqn{J} elements, where each element corresponds to a mixture component
#'   and contains a list of \eqn{T - 1} elements (one per time step).
#'   Each element in the time-step list is a list of \eqn{K} matrices (one for each past outcome category \eqn{y_{t-1}}),
#'   with dimensions \eqn{N \times (K-\mathrm{minimal})}, representing the score contributions for each transition \eqn{y_{t-1} \to y_t}.
#'
#' @examples
#' outcomes <- 1:2
#' T_max <- 3
#' J <- 2
#' model <- generate_transition_model(outcomes = outcomes, J = J, T_max = T_max)
#' sample <- generate_sample(model, N = 100)
#' Ps <- model$Ps_control  # List of transition matrices for each component
#' y_indices_matrix <- sample$y
#' posteriors <- sample$posteriors
#'
#' scores <- score_Ps(y_indices_matrix, Ps, posteriors)
#' print(length(scores))  # Should be J
#' print(length(scores[[1]]))  # Should be T-1 (time steps)
#' print(dim(scores[[1]][[1]][[1]]))  # Should be (N, K)
#'
#' @export
score_Ps <- function(y_indices_matrix, Ps, posteriors, minimal = TRUE) {
  # Check if every P_j in Ps have the same length
  if (!all(sapply(Ps, function(P_j) length(P_j) == length(Ps[[1]])))) {
    stop("All P_j in Ps must have the same length")
  }

  # Initialize a list to store score matrices
  scores <- list()

  # Iterate over all j
  for (j in 1:length(Ps)) {
    scores[[j]] <- score_P_j(y_indices_matrix, Ps[[j]], posteriors[,j], j_index = j,
                             minimal = minimal)
  }

  # Return the list of list of score matrices for all j
  scores
}

#' Compute Score Function for Mixture Component Priors
#'
#' This function calculates the score function for the prior probabilities of mixture components
#' in a mixture model, given sequences of observed outcomes and posterior probabilities.
#'
#' @param priors A numeric vector of length \eqn{J}, representing prior probabilities for each of \eqn{J} mixture components.
#' @param posteriors A numeric \eqn{N \times J} matrix, where each entry represents
#'   the posterior probability of observation \eqn{i} belonging to component \eqn{j}.
#' @param minimal A logical value. If TRUE, returns \eqn{N \times (J-1)} matrices.
#'
#' @return A numeric \eqn{N \times (J-\mathrm{minimal})} matrix, where each row corresponds to an observation, and each column
#'   represents the score function for a given mixture component prior.
#'
#' @examples
#' N <- 100
#' J <- 3
#' priors <- c(0.3, 0.5, 0.2)
#' posteriors <- matrix(runif(N * J), nrow = N, ncol = J)
#' posteriors <- posteriors / rowSums(posteriors)  # Normalize rows to sum to 1
#' y_indices_matrix <- matrix(sample(1:3, N * 4, replace = TRUE), nrow = N)
#'
#' scores <- score_priors(y_indices_matrix, priors, posteriors)
#' print(dim(scores))  # Should be (N, J)
#'
#' @export
score_priors <- function(priors, posteriors, minimal = TRUE) {
  # Check if length(priors) = ncol(posteriors)
  if (length(priors) != ncol(posteriors)) {
    stop("Length of priors must be equal to the number of columns in posteriors")
  }

  # Check if any of priors is zero
  if (any(priors == 0)) {
    stop("Prior probabilities must be non-zero")
  }

  # Compute the score for each prior
  J <- length(priors)
  if (J - minimal < 1) {
    # Return a matrix of zeros (no variation)
    return(matrix(0, nrow = nrow(posteriors), ncol = 0))
  }

  scores <- matrix(sapply(1:J, function(j) posteriors[,j] / priors[j] - posteriors[,J] / priors[J]),
                   ncol = J)

  # scores[,J] defined to be 1 - sum of all other scores
  if (J > 1) {
    scores[,J] <- 1 - rowSums(scores)
  }

  # Return the scores
  scores <- matrix(scores[,(1:(J-minimal))], ncol = (J - minimal))
  colnames(scores) <- paste0("prior_j", 1:(J-minimal))

  scores
}
