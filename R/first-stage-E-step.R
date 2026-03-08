#' Perform the E-Step of an EM Algorithm
#'
#' This function executes the E-Step in an EM algorithm, computing the
#' posterior probabilities of latent group membership for each unit
#' based on observed state sequences and current model parameters.
#'
#' @param transition_model A transition model with specified outcomes and number of latent groups.
#' @param data_for_est A list containing the state sequences for control and treated units,
#'              including `y_indices_control` and `y_indices_treated`.
#'
#' @export
E_step <- function(transition_model, data_for_est) {
  is_treated <- data_for_est$is_treated
  is_control <- data_for_est$is_control
  posteriors <- matrix(0, nrow = data_for_est$n, ncol = get_J(transition_model))

  # Initial P(Y_1 | D, Z) distribution
  p_y1d_splitted <- split_p_y1d(transition_model$p_y1d)

  # Extract parameters from transition model
  posteriors_control <- transitions_to_posteriors(data_for_est$y_indices_for_e_step_control,
                                                  p_y1d_splitted$p_y1d_control,
                                                  transition_model$Ps_control,
                                                  transition_model$priors)

  # Compute posteriors for the treated;
  # Note that transition matrices are the same for treated and control in pretreatment periods
  posteriors_treated <- transitions_to_posteriors(data_for_est$y_indices_for_e_step_treated,
                                                  p_y1d_splitted$p_y1d_treated,
                                                  transition_model$Ps_treated,
                                                  transition_model$priors)

  posteriors[is_control,] <- posteriors_control
  posteriors[is_treated,] <- posteriors_treated
  return(as.matrix(posteriors))
}

#' Split Joint Initial Distribution by Treatment Status
#'
#' Splits the joint pmf p(Y_1, D) into separate distributions for
#' control and treated units.
#'
#' @param p_y1d A list of matrices (one per latent group) where columns
#'   represent treatment status (1=control, 2=treated).
#' @return A list with `p_y1d_control` and `p_y1d_treated`, each a list
#'   of vectors (one per latent group).
#' @keywords internal
split_p_y1d <- function(p_y1d) {
  list(p_y1d_control = lapply(p_y1d, function(p_y1d_j) p_y1d_j[,1]),
       p_y1d_treated = lapply(p_y1d, function(p_y1d_j) p_y1d_j[,2]))
}

