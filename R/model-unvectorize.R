#' Unvectorize Transition Model Components
#'
#' This function reconstructs a transition model from its vectorized components.
#'
#' @param vectorized_model A list of numeric vectors representing the vectorized transition model components.
#' @param model_structure A list representing the original structure of the transition model (to infer correct dimensions).
#' @param minimal Logical; if TRUE, reconstructs missing probabilities in PMFs and transition matrices.
#' @return A list representing the original transition model.
#' @export
unvectorize_transition_model <- function(vectorized_model, model_structure, minimal = TRUE) {
  # Recover J
  priors <- unvectorize_pmf(vectorized_model$priors, minimal)
  J <- length(priors)

  # Recover Ps_nrows, from model_struture
  Ps_nrows <- sapply(model_structure$Ps_control_empirical[[1]], nrow)
  reconstructed_model <- list(
    priors = priors,
    priors_treated = unvectorize_pmf(vectorized_model$priors_treated, minimal),
    pmfs_initial_control = unvectorize_pmfs(vectorized_model$pmfs_initial_control,
                                            length(model_structure$pmfs_initial_control), minimal),
    pmfs_initial_treated = unvectorize_pmfs(vectorized_model$pmfs_initial_treated,
                                            length(model_structure$pmfs_initial_treated), minimal),
    pmfs_pre_treatment_control = unvectorize_pmfs(vectorized_model$pmfs_pre_treatment_control,
                                                  length(model_structure$pmfs_pre_treatment_control), minimal),
    pmfs_pre_treatment_treated = unvectorize_pmfs(vectorized_model$pmfs_pre_treatment_treated,
                                                  length(model_structure$pmfs_pre_treatment_treated), minimal),
    Ps_control_empirical = unvectorize_transition_matrices_list(vectorized_model$Ps_control_empirical,
                                                           nrows = Ps_nrows,
                                                           list_length = J,
                                                           minimal),
    Ps_treated_empirical = unvectorize_transition_matrices_list(vectorized_model$Ps_treated_empirical,
                                                           nrows = Ps_nrows,
                                                           list_length = J,
                                                           minimal)
  )

  return(reconstructed_model)
}
