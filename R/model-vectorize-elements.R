
#' Vectorize Probability Mass Function (PMF)
#'
#' This function vectorizes a probability mass function by optionally removing the last element.
#'
#' @param pmf A numeric vector representing the PMF.
#' @param minimal Logical; if TRUE, removes the last element.
#' @return A numeric vector representing the vectorized PMF.
#' @export
vectorize_pmf <- function(pmf, minimal = TRUE) {
  if (!minimal) {
    return(pmf)
  }
  if (length(pmf) < 2) {
    return(numeric(0))
  }
  return(pmf[1:(length(pmf)-1)])
}

#' Vectorize Multiple PMFs
#'
#' This function applies `vectorize_pmf` to a list of PMFs.
#'
#' @param pmfs A list of numeric vectors representing PMFs.
#' @param minimal Logical; if TRUE, applies minimal vectorization.
#' @return A list of vectorized PMFs.
#' @export
vectorize_pmfs <- function(pmfs, minimal = TRUE) {
  c(sapply(pmfs, vectorize_pmf, minimal = minimal))
}

#' Vectorize Transition Matrix
#'
#' This function vectorizes a transition matrix by applying `vectorize_pmf` to each row.
#'
#' @param transition_matrix A numeric matrix representing the transition matrix.
#' @param minimal Logical; if TRUE, applies minimal vectorization.
#' @return A numeric vector representing the vectorized transition matrix.
#' @export
vectorize_transition_matrix <- function(transition_matrix, minimal = TRUE) {
  c(apply(transition_matrix, 1, vectorize_pmf, minimal = minimal))
}

#' Vectorize Multiple Transition Matrices
#'
#' This function applies `vectorize_transition_matrix` to a list of transition matrices.
#'
#' @param transition_matrices A list of numeric matrices representing transition matrices.
#' @param minimal Logical; if TRUE, applies minimal vectorization.
#' @return A list of vectorized transition matrices.
#' @export
vectorize_transition_matrices <- function(transition_matrices, minimal = TRUE) {
  return(c(sapply(transition_matrices, vectorize_transition_matrix, minimal = minimal)))
}

#' Vectorize Multiple Transition Matrix Lists
#'
#' This function applies `vectorize_transition_matrix` to each transition matrix in a list of lists.
#' It flattens all transition matrices into a single vector.
#'
#' @param transition_matrices_list A list of lists, where each inner list contains transition matrices.
#' @param minimal Logical; if TRUE, applies minimal vectorization by removing the last column.
#' @return A numeric vector representing all vectorized transition matrices from the input list.
#' @export
vectorize_transition_matrices_list <- function(transition_matrices_list, minimal = TRUE) {
  if (length(transition_matrices_list) == 0) {
    return(list())
  }

  c(sapply(transition_matrices_list, vectorize_transition_matrices, minimal = minimal))
}

