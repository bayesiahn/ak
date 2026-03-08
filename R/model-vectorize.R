#' Vectorize Transition Model Components
#'
#' This function extracts and vectorizes key elements from a transition model,
#' returning them as a list of numeric vectors.
#'
#' @param model A list representing the transition model, containing elements such as
#'        priors, priors_treated, probability mass functions (pmfs), and transition matrices.
#' @param minimal Logical; if TRUE, applies minimal vectorization by removing the last column
#'        of transition matrices and omitting the last probability in PMFs.
#'
#' @return A list of numeric vectors containing the vectorized elements of the transition model.
#' @export
vectorize_transition_model <- function(model, minimal = TRUE) {
  # Extract elements and convert them to vectors
  vectorized_elements <- list(
    priors = vectorize_pmf(model$priors, minimal),
    priors_treated = vectorize_pmf(model$priors_treated, minimal),
    pmfs_initial_control = vectorize_pmfs(model$pmfs_initial_control, minimal),
    pmfs_initial_treated = vectorize_pmfs(model$pmfs_initial_treated, minimal),
    pmfs_pre_treatment_control = vectorize_pmfs(model$pmfs_pre_treatment_control, minimal),
    pmfs_pre_treatment_treated = vectorize_pmfs(model$pmfs_pre_treatment_treated, minimal),
    Ps_control_empirical = vectorize_transition_matrices_list(model$Ps_control_empirical, minimal),
    Ps_treated_empirical = vectorize_transition_matrices_list(model$Ps_treated_empirical, minimal)
  )

  return(vectorized_elements)
}

#' Vectorize Transition Model Components
#'
#' This function extracts and vectorizes key elements from a transition model,
#' returning them as a vector
#'
#' @param model A list representing the transition model, containing elements such as
#'        priors, priors_treated, probability mass functions (pmfs), and transition matrices.
#' @param minimal Logical; if TRUE, applies minimal vectorization by removing the last column
#'        of transition matrices and omitting the last probability in PMFs.
#'
#' @return A numeric vector consisting of vectorized elements of the transition model.
#' @export
completely_vectorize_transition_model <- function(model, minimal = TRUE) {
  unlist(vectorize_transition_model(model, minimal))
}
