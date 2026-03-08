#' Create Block Diagonal Matrix from Two Matrices
#'
#' Combines two matrices into a block diagonal matrix. This is a simple
#' implementation to avoid dependency on the magic package.
#'
#' @param A First matrix (top-left block).
#' @param B Second matrix (bottom-right block).
#' @return A block diagonal matrix with A in the top-left and B in the bottom-right.
#' @keywords internal
block_diagonal <- function(A, B) {
  rbind(
    cbind(A, matrix(0, nrow = nrow(A), ncol = ncol(B))),
    cbind(matrix(0, nrow = nrow(B), ncol = ncol(A)), B)
  )
}

#' Combine List of Matrices Column-wise
#'
#' This function takes a list of matrices or a list of lists of matrices and binds them column-wise.
#' If the input is a list of matrices, it simply applies `cbind`. If the input is a list of lists of matrices,
#' it recursively applies the function to each sublist and then binds them column-wise.
#'
#' @param l A list of matrices or a list of lists of matrices.
#' @return A combined matrix with elements bound column-wise.
#' @keywords internal
#' @examples
#' mat_list <- list(matrix(1:4, 2, 2), matrix(5:8, 2, 2))
#' cbind_list_of_matrices(mat_list)
cbind_list_of_matrices <- function(l) {
  if (all(sapply(l, function(el) is.matrix(el) || is.numeric(el)))) {
    do.call(cbind, l)
  } else if (all(sapply(l, is.list))) {
    do.call(cbind, lapply(l, cbind_list_of_matrices))
  } else {
    stop("Input must be a list of matrices or a list of lists of matrices.")
  }
}

#' Convert Model to Score Matrix
#'
#' This function computes the score matrix from a given model and data. It extracts relevant probability mass functions
#' and transition matrices from the model, vectorizes them, and combines them into a single score matrix.
#'
#' @param y Response variable.
#' @param g Grouping variable.
#' @param model A fitted model containing priors and transition matrices.
#' @return A score matrix representing the structured parameter estimates.
#' @keywords internal
model_to_score_matrix  <- function(y, g, model) {
  scores <- model_to_scores(y, g, model)
  J <- length(model$priors)

  vectorized_elements_common <- list()
  vectorized_elements_treated <- list(
    pmfs_initial_treated = cbind_list_of_matrices(scores$pmfs_initial_treated),
    Ps_treated_empirical = cbind_list_of_matrices(scores$Ps_treated_empirical)
  )
  vectorized_elements_control <- list(
    pmfs_initial_control = cbind_list_of_matrices(scores$pmfs_initial_control),
    Ps_control_empirical = cbind_list_of_matrices(scores$Ps_control_empirical)
  )

  if (J > 1) {
    vectorized_elements_common <- list(priors = scores$priors)
    vectorized_elements_treated <- c(vectorized_elements_treated, list(priors_treated = scores$priors_treated))
  }

  scores_control <- do.call(cbind, vectorized_elements_control)
  scores_treated <- do.call(cbind, vectorized_elements_treated)
  colnames(scores_control) <- paste0("d0_", colnames(scores_control))
  colnames(scores_treated) <- paste0("d1_", colnames(scores_treated))

  scores_diagonal <- block_diagonal(scores_control, scores_treated)

  cbind(scores_diagonal, do.call(cbind, vectorized_elements_common))
}

#' Convert Model to Vector for Confidence Interval Computation
#'
#' This function transforms the model parameters into a single vector, suitable for confidence interval estimation.
#' It extracts and flattens probability mass functions, transition matrices, and priors.
#'
#' @param model A fitted model.
#' @return A numeric vector of model parameters.
#' @keywords internal
model_to_vector_for_ci <- function(model) {
  J <- length(model$priors)

  vectorized_model_for_ci <- c(
    vectorize_pmfs(model$pmfs_initial_control),
    vectorize_transition_matrices_list(model$Ps_control_empirical),
    vectorize_pmfs(model$pmfs_initial_treated),
    vectorize_transition_matrices_list(model$Ps_treated_empirical)
  )

  if (J > 1) {
    vectorized_model_for_ci <- c(
      vectorized_model_for_ci,
      model$priors_treated[1:(J-1)],
      model$priors[1:(J-1)]
    )
  }

  vectorized_model_for_ci
}

#' Compute Asymptotic Variance from Score Matrix
#'
#' This function estimates the asymptotic variance of the parameter estimates based on the score matrix.
#'
#' @param score_matrix A numeric score matrix.
#' @return The asymptotic variance matrix.
#' @keywords internal
score_matrix_to_variance_asymptotic <- function(score_matrix) {
  N <- nrow(score_matrix)
  solve(t(score_matrix) %*% score_matrix / N)
}

#' Compute Asymptotic Standard Errors from Score Matrix
#'
#' This function derives standard errors from the asymptotic variance of the score matrix.
#'
#' @param score_matrix A numeric score matrix.
#' @return A vector of standard errors for each parameter.
#' @keywords internal
score_matrix_to_se_asymptotic <- function(score_matrix) {
  V <- score_matrix_to_variance_asymptotic(score_matrix)
  N <- nrow(score_matrix)
  sqrt(diag(V) / N)
}

#' Compute Asymptotic Confidence Intervals for Model Parameters
#'
#' This function computes confidence intervals for the model parameters using an asymptotic approximation.
#' It constructs a score matrix, computes standard errors, and applies the normal quantile function to generate bounds.
#'
#' @param model A fitted model.
#' @param y Outcome data (N-length vector)
#' @param g Treatment cohort data (N-length vector)
#' @param c Cluster indicator data (N-length vector, optional, default: NULL).
#' @param lower_percentile Lower percentile for confidence interval (default: 0.025).
#' @param upper_percentile Upper percentile for confidence interval (default: 0.975).
#' @return A list containing lower and upper confidence intervals.
#' @keywords internal
model_to_cis_asymptotic <- function(model, y, g, c = NULL, lower_percentile = 0.025, upper_percentile = 0.975) {
  score_matrix <- model_to_score_matrix(y, g, model)
  se <- score_matrix_to_se_asymptotic(score_matrix)
  z <- qnorm(c(lower_percentile, upper_percentile))
  model_vector <- model_to_vector_for_ci(model)

  N <- nrow(score_matrix)
  ci_lower <- model_vector + z[1] * se
  ci_upper <- model_vector + z[2] * se

  list(lower = c(ci_lower), upper = c(ci_upper))
}
