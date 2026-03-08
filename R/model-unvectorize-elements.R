#' Unvectorize Probability Mass Function (PMF)
#'
#' This function reconstructs a probability mass function from a vectorized form.
#' If `minimal = TRUE`, it ensures the last probability is reconstructed so the sum equals 1.
#'
#' @param pmf_vector A numeric vector representing a vectorized PMF.
#' @param minimal Logical; if TRUE, reconstructs the last probability to ensure sum = 1.
#' @return A numeric vector representing the original PMF.
#' @export
unvectorize_pmf <- function(pmf_vector, minimal = TRUE) {
  if (!minimal) {
    return(pmf_vector)
  }

  # If empty, return a single 1.0 probability
  if (length(pmf_vector) == 0) {
    return(c(1.0))
  }

  # Compute the last probability to ensure the sum equals 1
  last_prob <- 1 - sum(pmf_vector)
  return(c(pmf_vector, last_prob))
}

#' Unvectorize Transition Matrix
#'
#' This function reconstructs a transition matrix from its vectorized form.
#' If `minimal = TRUE`, it ensures that the last column is computed so row sums equal 1.
#'
#' @param vectorized_matrix A numeric vector representing a row-wise vectorized transition matrix.
#' @param nrow The number of rows in the original transition matrix, i.e.,
#' the number of current states (vs future)
#' @param minimal Logical; if TRUE, reconstructs the last column to ensure row sums equal 1.
#' @return A numeric matrix representing the original transition matrix.
#' @export
unvectorize_transition_matrix <- function(vectorized_matrix, nrow, minimal = TRUE) {
  if (!minimal) {
    return(matrix(vectorized_matrix, nrow = nrow, byrow = TRUE))
  }

  # Determine number of columns (last column is missing)
  n_cols <- (length(vectorized_matrix) / nrow) + 1

  # Reshape into matrix with missing last column
  restored_matrix <- matrix(vectorized_matrix, nrow = nrow, byrow = TRUE)

  # Compute last column to ensure row sums equal 1 if minimal
  if (minimal) {
    last_col <- 1 - rowSums(restored_matrix)
  }

  # Reconstruct full transition matrix
  full_matrix <- cbind(restored_matrix, last_col)
  dimnames(full_matrix) <- NULL
  return(full_matrix)
}

#' Unvectorize Transition Matrix
#'
#' This function reconstructs a transition matrix from its vectorized form.
#' If `minimal = TRUE`, it ensures that the last column is computed so row sums equal 1.
#'
#' @param vectorized_matrix A numeric vector representing a row-wise vectorized transition matrix.
#' @param nrow The number of rows in the original transition matrix, i.e.,
#' the number of current states (vs future)
#' @param minimal Logical; if TRUE, reconstructs the last column to ensure row sums equal 1.
#' @return A numeric matrix representing the original transition matrix.
#' @export
unvectorize_transition_matrices <- function(vectorized_matrices, nrows, minimal = TRUE) {
  if (!minimal) {
    return(matrix(vectorized_matrix, nrow = nrow, byrow = TRUE))
  }

  # Determine number of columns (last column is missing)
  n_cols <- (length(vectorized_matrix) / nrow) + 1

  # Reshape into matrix with missing last column
  restored_matrix <- matrix(vectorized_matrix, nrow = nrow, byrow = TRUE)

  # Compute last column to ensure row sums equal 1 if minimal
  if (minimal) {
    last_col <- 1 - rowSums(restored_matrix)
  }

  # Reconstruct full transition matrix
  full_matrix <- cbind(restored_matrix, last_col)
  dimnames(full_matrix) <- NULL
  return(full_matrix)
}

#' Unvectorize Multiple Transition Matrices
#'
#' This function reconstructs a list of transition matrices from a single vector.
#' If `minimal = TRUE`, it ensures the last column of each matrix is reconstructed so row sums equal 1.
#'
#' @param vectorized_matrices A numeric vector representing multiple row-wise vectorized transition matrices.
#' @param nrows A numeric vector indicating the number of rows in each original transition matrix.
#' @param minimal Logical; if TRUE, reconstructs the last column to ensure row sums equal 1.
#' @return A list of numeric matrices representing the original transition matrices.
#' @export
unvectorize_transition_matrices <- function(vectorized_matrices, nrows, minimal = TRUE) {
  if (length(nrows) == 0) {
    return(list())
  }

  matrices <- list()
  start_idx <- 1

  for (i in seq_along(nrows)) {
    ncols_i <- ifelse(i < length(nrows), nrows[(i+1)] - minimal,
                      (length(vectorized_matrices) - start_idx + 1) / nrows[i])
    i_elements_length <- nrows[i] * ncols_i
    end_idx <- start_idx + i_elements_length - 1

    matrix_vector <- vectorized_matrices[start_idx:end_idx]
    matrices[[i]] <- unvectorize_transition_matrix(matrix_vector, nrow = nrows[i], minimal = minimal)

    start_idx <- end_idx + 1
  }

  return(matrices)
}

#' Unvectorize Multiple PMF Lists
#'
#' This function reconstructs a list of probability mass function (PMF) lists from a single vector.
#'
#' @param pmf_list A numeric vector representing multiple concatenated vectorized PMFs.
#' @param list_length The number of lists to split the PMFs into.
#' @param minimal Logical; if TRUE, reconstructs the last probability to ensure sum = 1.
#' @return A list of lists, where each inner list contains the original PMFs.
#' @export
unvectorize_pmfs <- function(vectorized_pmfs, list_length = 1, minimal = TRUE) {
  if (length(vectorized_pmfs) == 0) {
    return(replicate(list_length, list()))
  }

  # Check if vectorized_pmfs is dividable by list_length
  if (length(vectorized_pmfs) %% list_length != 0) {
    stop("Length of vectorized_pmfs must be divisible by list_length")
  }

  # Split the vectorized PMFs into equal parts
  pmf_splits <- split(vectorized_pmfs,
                      ceiling(seq_along(vectorized_pmfs) / (length(vectorized_pmfs) / list_length)))

  # Apply unvectorize_pmf to each split
  return(unname(lapply(pmf_splits, unvectorize_pmf, minimal = minimal)))
}

#' Unvectorize Multiple Transition Matrix Lists
#'
#' This function reconstructs a list of transition matrix lists from a single vector.
#'
#' @param vectorized_matrices_list A numeric vector representing multiple concatenated vectorized transition matrices.
#' @param nrows A numeric vector indicating the number of rows in each original transition matrix.
#' @param list_length The number of lists to split the transition matrices into.
#' @param minimal Logical; if TRUE, reconstructs the last column to ensure row sums equal 1.
#' @return A list of lists, where each inner list contains the original transition matrices.
#' @export
unvectorize_transition_matrices_list <- function(vectorized_matrices_list, nrows,
                                                 list_length = 1, minimal = TRUE) {
  if (length(vectorized_matrices_list) == 0) {
    return(replicate(list_length, list()))
  }

  # Check if vectorized_matrices_list is dividable by list_length
  if (length(vectorized_matrices_list) %% list_length != 0) {
    stop("Length of vectorized_matrices_list must be divisible by list_length")
  }

  # Split the vectorized transition matrices into equal parts
  split_size <- length(vectorized_matrices_list) / list_length
  matrix_splits <- split(vectorized_matrices_list, ceiling(seq_along(vectorized_matrices_list) / split_size))

  # Apply unvectorize_transition_matrices to each split
  return(unname(lapply(matrix_splits,
    unvectorize_transition_matrices, nrows = nrows, minimal = minimal)))
}

