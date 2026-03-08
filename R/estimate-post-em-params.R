#' Estimate Post-Treatment Parameters
#'
#' This function estimates post-treatment parameters for the second stage of the
#' algorithm, including transition matrices for treated and control groups and
#' probability mass functions (PMFs) for states just before treatment.
#'
#' @param posteriors A matrix of posterior probabilities where each column
#'        corresponds to a latent group.
#' @param data_for_est A list containing processed data for estimation,
#'        including state sequences and indicators.
#' @param unit_weights (Optional) A numeric vector of weights for each unit.
#' @param opts (Optional) A list of options for the algorithm.
#'
#' @return A list containing estimated post-treatment parameters:
#'         PMFs for states just before treatment, and transition matrices for
#'         treated and control groups.
#' @export
estimate_post_em_params <- function(posteriors, data_for_est,
                                           unit_weights = NULL,
                                           opts = list()) {
  # Compute weights from posteriors; if hard classification is requested, convert posteriors to hard labels
  weights <- posteriors
  if (check_if_hard_classification(opts)) {
    weights <- hard_classification(posteriors)
  }

  # Validate and handle unit weights
  unit_weights <- check_unit_weights_validity(data_for_est$y_indices_matrix, unit_weights)

  # Extract the model parameters
  J <- ncol(posteriors)
  is_control <- data_for_est$is_control
  is_treated <- data_for_est$is_treated
  unit_weights_treated <- unit_weights[is_treated]
  unit_weights_control <- unit_weights[is_control]
  weights_control <- matrix(weights[is_control, ], ncol = J)
  weights_treated <- matrix(weights[is_treated, ], ncol = J)

  # Compute pmfs_pre_treatment_treated
  pmfs_pre_treatment_treated <- M_step_initial_outcome_dists(weights_treated,
    data_for_est$indicators_treated_just_before_treatment,
    unit_weights = unit_weights_treated)
  pmfs_pre_treatment_control <- M_step_initial_outcome_dists(weights_control,
    data_for_est$indicators_control_just_before_treatment,
    unit_weights = unit_weights_control)

  # Compute empirical Ps_control and Ps_treated (Markovian one-step transitions)
  Ps_control_empirical <- M_step_Ps(weights_control, data_for_est$transition_indicators_control,
                                    unit_weights = unit_weights_control)
  Ps_treated_empirical <- M_step_Ps(weights_treated, data_for_est$transition_indicators_treated,
                                    unit_weights = unit_weights_treated)

  # Compute non-Markovian transitions from the initial period (t=0) to each subsequent period
  # Note: Ps_*_from_0[[1]] should match Ps_*_empirical[[1]] by construction
  Ps_control_from_0 <- M_step_Ps(weights_control, data_for_est$transition_indicators_control_baseline,
                                 unit_weights = unit_weights_control)
  Ps_treated_from_0 <- M_step_Ps(weights_treated, data_for_est$transition_indicators_treated_baseline,
                                 unit_weights = unit_weights_treated)

  # Compute probability of treatment assignment and priors for treated groups
  treated_period <- max(data_for_est$g)
  prob_treated <- compute_prob_treated(data_for_est$g, unit_weights = unit_weights)
  priors_treated <- compute_priors_treated(weights, data_for_est$g, unit_weights = unit_weights)


  return(list(pmfs_pre_treatment_control = pmfs_pre_treatment_control,
              pmfs_pre_treatment_treated = pmfs_pre_treatment_treated,
              Ps_treated_empirical = Ps_treated_empirical, Ps_control_empirical = Ps_control_empirical,
              Ps_control_from_0 = Ps_control_from_0, Ps_treated_from_0 = Ps_treated_from_0,
              treated_period = treated_period,
              prob_treated = prob_treated, priors_treated = priors_treated,
              weights = weights, posteriors = posteriors,
              unit_weights = unit_weights))
}

#' Compute Probability of Treatment Assignment
#'
#' This function calculates the probability of treatment assignment by
#' computing the weighted mean of the treatment indicator.
#'
#' @param g A vector indicating the treatment status for each unit.
#' @param unit_weights A vector of individual weights for each unit. If NULL, all
#'  units are assumed to have equal weight.
#'
#' @return A numeric value representing the probability of treatment assignment.
#'
#' @examples
#' g <- sample(0:1, 10, replace = TRUE)
#' prob_treated <- compute_prob_treated(g)
#'
#' unit_weights <- runif(10)
#' prob_treated_weighted <- compute_prob_treated(g, unit_weights)
#'
#' @export
compute_prob_treated <- function(g, unit_weights = NULL) {
  # Validate and handle unit weights
  if (is.null(unit_weights)) {
    unit_weights <- rep(1, length(g))
  }
  if (length(g) == 0) {
    stop("Input vector 'g' is empty.")
  }
  if (length(g) != length(unit_weights)) {
    stop("Length of 'g' and 'unit_weights' must be the same.")
  }


  # Compute weighted mean
  unit_weights_treated <- unit_weights[g>0]
  sum(unit_weights_treated) / sum(unit_weights)
}

#' Compute Priors for Treated Groups
#'
#' This function calculates the mean posterior probabilities for the treated
#' units, effectively estimating the priors for the treated groups.
#'
#' @param weights An N by J matrix of individual weights of group membership.
#' @param g A vector indicating the treatment status for each unit.
#' @param unit_weights A vector of individual weights for each unit. If NULL, all
#'  units are assumed to have equal weight.
#'
#' @return A vector of mean posterior probabilities for each treated group.
#'
#' @examples
#' weights <- matrix(runif(20), nrow = 10, ncol = 2)
#' g <- sample(0:1, 10, replace = TRUE)
#' priors_treated <- compute_priors_treated(weights, g)
#'
#' @export
compute_priors_treated <- function(weights, g, unit_weights = NULL) {
  # Extract treated weights
  J <- ncol(weights)
  weights_treated <- matrix(weights[g > 0,], ncol = J)
  unit_weights_treated <- unit_weights[g > 0]

  # Compute priors for treated groups
  M_step_priors(weights_treated, unit_weights_treated)
}

#' Aggregate Average Treatment Effects (LTATTs) Across Groups
#'
#' This function aggregates LTATTs across different transition patterns, weighted by
#' their respective treated priors.
#'
#' @param ATT_list A list of ATT vectors for each transition pattern.
#' @param priors_treated A vector of treated priors for each transition pattern.
#'
#' @return A single vector of aggregated LTATTs.
#'
#' @examples
#' ATT_list <- list(c(0.1, 0.2), c(0.3, 0.4))
#' priors_treated <- c(0.5, 0.5)
#' aggregated_ATTs <- aggregate_ATT_by_latent_type(ATT_list, priors_treated)
#'
#' @export
aggregate_ATT_by_latent_type <- function(ATT_list, priors_treated) {
  # Check if every element of ATT_list has the same length
  if (length(unique(sapply(ATT_list, length))) > 1) {
    stop("LTATTs must have the same length.")
  }

  # Multiply every jth member of ATT_list by jth member priors_treated and sum
  J <- length(ATT_list)
  ATT <- sapply(1:length(ATT_list[[1]]), function(t) {
    sum(sapply(1:length(ATT_list), function(j) {
      ATT_list[[j]][t] * priors_treated[j]
    }))
  })

  ATT
}
