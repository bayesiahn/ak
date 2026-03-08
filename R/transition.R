#' Obtain the index of each element in y_t within Y_t$values
#'
#' @param y_t Numeric vector representing the outcomes for each unit.
#' @param Y_t List containing 'values', a sorted numeric vector, and 'type',
#' a character string indicating if the outcome space is 'discrete' or 'continuous'.
#' For continuous measures use quantile indices; first index such that
#' y_{it} <= Y_t$values is TRUE
#' @return Returns a numeric vector of indices indicating the position of
#' each element of y_t within Y_t$values.
#'
#' @examples
#' y_t <- c(1, 2, 3, 3)
#' Y_t <- list(values = c(1, 2, 3), type = 'discrete')
#' get_y_t_by_Y_t_index(y_t, Y_t) # should return c(1, 2, 3, 3)
#'
#' y_t <- c(0, 0.2, 0.4, 0.5)
#' Y_t <- list(values = c(0.1, 0.3, 0.5), type = 'continuous')
#' get_y_t_by_Y_t_index(y_t, Y_t) # should return c(1, 2, 3, 3)
#' @export
get_y_t_by_Y_t_index <- function(y_t, Y_t) {
  if (Y_t$type == "discrete") {
    return(match(y_t, Y_t$values))
  } else if (Y_t$type == "continuous") {
    indices <- findInterval(y_t, Y_t$values) + 1
    exact_matches <- y_t %in% Y_t$values
    indices[exact_matches] <- match(y_t[exact_matches], Y_t$values)
    return(indices)
  } else {
    stop("Unknown type specified for Y_t")
  }
}

#' Compute Transition Matrix Between Two Outcome Spaces
#'
#' This function computes the transition matrix that represents the
#' probability of transitioning from an outcome in `Y_t_now` to
#' an outcome in `Y_t_forward`.
#'
#' @param y_matrix A matrix of size N by 2 representing the outcome histories.
#'        The ith row represents the outcome history of unit i for now and forward, respectively.
#' @param Y_t_now A list containing outcome space details for the 'now' period.
#'        It should contain `values`, a numeric vector, and `type`, a character string.
#' @param Y_t_forward A list containing outcome space details for the 'forward' period.
#'        It should contain `values`, a numeric vector, and `type`, a character string.
#' @return A matrix where the (i,j)th element represents the probability of
#'         transitioning from the ith element of `Y_t_now$values` to the jth element
#'         of `Y_t_forward$values`.
#' @examples
#' y_matrix <- matrix(c(1, 2, 3, 2, 1, 3), ncol = 2)
#' Y_t_now <- list(values = c(1, 2, 3), type = 'discrete')
#' Y_t_forward <- list(values = c(1, 2, 3), type = 'discrete')
#' transition_matrix <- get_transition_matrix(y_matrix, Y_t_now, Y_t_forward)
#' print(transition_matrix)
#' @export
get_transition_matrix <- function(y_matrix, Y_t_now, Y_t_forward) {
  N <- length(y_matrix[,1])
  y_t_now <- y_matrix[,1]
  y_t_forward <- y_matrix[,2]

  Y_t_now_indices <- 1:length(Y_t_now$names)
  Y_t_forward_indices <- 1:length(Y_t_forward$names)

  y_t_now_by_Y_t_index <- match(y_t_now, Y_t_now$names)
  y_t_forward_by_Y_t_index <- match(y_t_forward, Y_t_forward$names)

  transition_matrix <- matrix(nrow = length(Y_t_now_indices),
                              ncol = length(Y_t_forward_indices))
  for (l in Y_t_now_indices) {
    denom <- sum(y_t_now_by_Y_t_index == l)
    for (k in Y_t_forward_indices) {
      transition_matrix[l,k] = sum(
        (y_t_forward_by_Y_t_index == k) * (y_t_now_by_Y_t_index == l)) /
        denom
    }
  }

  rownames(transition_matrix) <- get_Y_t_names(Y_t_now)
  colnames(transition_matrix) <- get_Y_t_names(Y_t_forward)

  return(transition_matrix)
}

#' Apply Transition Matrices to a PMF
#'
#' This function computes the resulting PMF after sequentially applying a list of transition matrices
#' to an initial PMF.
#'
#' @param initial_pmf A numeric vector representing the initial PMF.
#' @param transition_matrices A list of matrices, where each matrix represents a transition probability matrix.
#'
#' @return A numeric vector representing the resulting PMF after applying all the transition matrices.
#' @examples
#' # Example initial PMF
#' initial_pmf <- c(0.5, 0.3, 0.2)
#'
#' # Example transition matrices
#' transition_matrices <- list(
#'   matrix(c(0.8, 0.1, 0.1, 0.2, 0.6, 0.2, 0.3, 0.3, 0.4), nrow = 3, byrow = TRUE),
#'   matrix(c(0.7, 0.2, 0.1, 0.2, 0.5, 0.3, 0.4, 0.4, 0.2), nrow = 3, byrow = TRUE)
#' )
#'
#' # Compute the resulting PMF
#' resulting_pmf <- apply_transition_matrices(initial_pmf, transition_matrices)
#' print(resulting_pmf)
#'
#' @export
apply_transition_matrices <- function(initial_pmf, transition_matrices) {
  # Validate inputs
  if (!is.numeric(initial_pmf) || !all(initial_pmf >= 0)) {
    stop("initial_pmf must be a numeric vector with non-negative values.")
  }
  if (!is.list(transition_matrices) || !all(sapply(transition_matrices, is.matrix))) {
    stop("transition_matrices must be a list of matrices.")
  }

  # Ensure the PMF sums to 1
  initial_pmf <- t(initial_pmf / sum(initial_pmf))

  # Sequentially apply each transition matrix
  resulting_pmf <- initial_pmf
  for (matrix in transition_matrices) {
    if (ncol(matrix) != length(resulting_pmf)) {
      stop("Each transition matrix must have a number of columns equal to the length of the PMF.")
    }
    resulting_pmf <- resulting_pmf %*% matrix
  }

  # Ensure the result is normalized
  resulting_pmf <- as.numeric(resulting_pmf)
  resulting_pmf / sum(resulting_pmf)
}


#' #' Compute Transition Matrix Between Two Outcome Spaces
#' #'
#' #' This function computes the transition matrix that represents the
#' #' probability of transitioning from an outcome in `Y_t_now` to
#' #' an outcome in `Y_t_forward`.
#' #'
#' #' @param y_matrix A matrix of size N by 2 representing the outcome histories.
#' #'        The ith row represents the outcome history of unit i for now and forward, respectively.
#' #' @return A matrix where the element [i, j] represents the probability of
#' #'         transitioning from the ith element of `Y_t_now$values` to the jth element
#' #'         of `Y_t_forward$values`.
#' #' @examples
#' #' y_matrix <- matrix(c(1, 2, 3, 2, 1, 3), ncol = 2)
#' #' transition_matrix <- get_transition_matrix(y_matrix)
#' #' print(transition_matrix)
#' #' @export
#' get_transition_matrix <- function(y_matrix) {
#'   get_transition_matrix(y_matrix, generate_Y_t(y_matrix[,1]), generate_Y_t(y_matrix[,2]))
#' }
