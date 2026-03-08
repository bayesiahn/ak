#' Compute Score Functions from a Model
#'
#' This function computes a set of score functions for a given model,
#' including priors, transition probabilities, and probability mass functions (PMFs)
#' for both treated and control groups.
#'
#' @param y A matrix of observed outcomes, where each row represents an individual and
#'   each column represents a time period.
#' @param g A numeric vector indicating treatment assignment, where values greater than
#'   zero indicate treated units and values of zero indicate control units.
#' @param model A list containing the following model components:
#'   - `priors`: Prior probabilities for mixture components.
#'   - `priors_treated`: Prior probabilities for mixture components in the treated group.
#'   - `Ps_control`: List of transition matrices for the control group.
#'   - `Ps_treated`: List of transition matrices for the treated group.
#'   - `Ps_control_empirical`: Empirical transition matrices for the control group.
#'   - `Ps_treated_empirical`: Empirical transition matrices for the treated group.
#'   - `pmfs_initial_control`: Initial probability mass functions for the control group.
#'   - `pmfs_initial_treated`: Initial probability mass functions for the treated group.
#'   - `pmfs_pre_treatment_control`: Pre-treatment PMFs for the control group.
#'   - `pmfs_pre_treatment_treated`: Pre-treatment PMFs for the treated group.
#'   - `treated_period`: Time period when treatment starts.
#'   - `Y`: List containing outcome names.
#'
#' @return A list of computed scores with the following elements:
#'   - `g`: Treatment assignment vector.
#'   - `Y`: Outcome names from the model.
#'   - `treated_period`: The period when treatment starts.
#'   - `priors`: Scores for mixture component priors.
#'   - `priors_treated`: Scores for treated priors.
#'   - `Ps_control`: Scores for transition probabilities in the control group.
#'   - `Ps_treated`: Scores for transition probabilities in the treated group.
#'   - `Ps_control_empirical`: Scores for empirical transition probabilities in control group.
#'   - `Ps_treated_empirical`: Scores for empirical transition probabilities in treated group.
#'   - `pmfs_initial_control`: Scores for initial PMFs in control group.
#'   - `pmfs_initial_treated`: Scores for initial PMFs in treated group.
#'   - `pmfs_pre_treatment_control`: Scores for pre-treatment PMFs in control group.
#'   - `pmfs_pre_treatment_treated`: Scores for pre-treatment PMFs in treated group.
#'
#' @examples
#' set.seed(123)
#' N <- 100
#' T_max <- 4
#' J <- 2
#' g <- sample(c(0, 1), N, replace = TRUE)  # Treatment assignment
#' y <- matrix(sample(1:3, N * T_max, replace = TRUE), nrow = N, ncol = T_max)
#' model <- generate_transition_model(outcomes = 1:3, J = J, T_max = T_max)
#'
#' scores <- model_to_scores(y, g, model)
#' print(names(scores))  # List of score components
#'
#' @export
model_to_scores <- function(y, g, model) {
  # Extract data
  treated_indices <- which(g > 0)
  treated_period <- model$treated_period
  T_max <- ncol(y)
  y_indices <- y_to_y_indices(y, model$Y$names)
  y_indices_control <- matrix(y_indices[-treated_indices,], nrow = sum(!(g > 0)))
  y_indices_treated <- matrix(y_indices[treated_indices,], nrow = sum(g > 0))

  # Compute posteriors for the control
  posteriors_control <- transitions_to_posteriors(y_indices_control,
                                                  model$pmfs_initial_control,
                                                  model$Ps_control,
                                                  model$priors)

  # Compute posteriors for the treated;
  # Note that transition matrices are the same for treated and control in pretreatment periods
  posteriors_treated <- transitions_to_posteriors(y_indices_treated,
                                                  model$pmfs_initial_treated,
                                                  model$Ps_treated,
                                                  model$priors)
  # Combine posteriors
  posteriors <- matrix(NA, nrow=length(g), ncol=ncol(posteriors_control))
  posteriors[-treated_indices,] <- posteriors_control
  posteriors[treated_indices,] <- posteriors_treated

  # Compute scores
  scores <- list()
  scores$g <- g
  scores$model <- model
  scores$Y <- model$Y
  scores$prob_treated <- model$prob_treated
  scores$treated_period <- model$treated_period
  scores$priors <- score_priors(model$priors, posteriors)
  scores$priors_treated <- score_priors(model$priors_treated, posteriors_treated)
  scores$pmfs_initial_control <- score_pmfs(y_indices_control[,1], model$pmfs_initial_control, posteriors_control, 1)
  scores$pmfs_initial_treated <- score_pmfs(y_indices_treated[,1], model$pmfs_initial_treated, posteriors_treated, 1)
  scores$pmfs_pre_treatment_control <- score_pmfs(y_indices_control[,(treated_period-1)], model$pmfs_pre_treatment_control, posteriors_control, treated_period - 1)
  scores$pmfs_pre_treatment_treated <- score_pmfs(y_indices_treated[,(treated_period-1)], model$pmfs_pre_treatment_treated, posteriors_treated, treated_period - 1)
  # Ps_control=Ps_treated is the same for treated and control in pre-treatment periods
  scores$Ps_treated_empirical <- score_Ps(y_indices_treated, model$Ps_treated_empirical, posteriors_treated)
  scores$Ps_control_empirical <- score_Ps(y_indices_control, model$Ps_control_empirical, posteriors_control)
  scores$Ps_treated <- score_Ps(y_indices, model$Ps_treated, posteriors)
  scores$Ps_control <- score_Ps(y_indices, model$Ps_control, posteriors)

  # Return scores
  scores
}

#' Convert Score List to Transition Model
#'
#' This function converts a list of score functions into a structured transition model,
#' which includes estimated priors, transition probabilities, and probability mass functions (PMFs).
#'
#' @param scores A list containing score matrices for priors, transition probabilities,
#'   and probability mass functions. The expected structure of `scores` includes:
#'   - `g`: Treatment assignment vector.
#'   - `priors`: A numeric score vector for prior probabilities.
#'   - `priors_treated`: A numeric score vector for prior probabilities in the treated group.
#'   - `Ps_control`: A list of transition score matrices for the control group.
#'   - `Ps_treated`: A list of transition score matrices for the treated group.
#'   - `Ps_control_empirical`: A list of empirical transition score matrices for the control group.
#'   - `Ps_treated_empirical`: A list of empirical transition score matrices for the treated group.
#'   - `pmfs_initial_control`: A list of score matrices for initial PMFs in the control group.
#'   - `pmfs_initial_treated`: A list of score matrices for initial PMFs in the treated group.
#'   - `pmfs_pre_treatment_control`: A list of score matrices for pre-treatment PMFs in the control group.
#'   - `pmfs_pre_treatment_treated`: A list of score matrices for pre-treatment PMFs in the treated group.
#'   - `treated_period`: The period when treatment starts.
#'   - `Y`: A list containing outcome names.
#' @param weights An optional numeric vector of weights for each observation, which is applied
#'   when computing priors, transition probabilities, and PMFs.
#'
#' @return A list representing a structured transition model with the following elements:
#'   - `priors`: Estimated prior probabilities.
#'   - `priors_treated`: Estimated prior probabilities for the treated group.
#'   - `prob_treated`: Estimated probability of treatment.
#'   - `Ps_control`: Estimated transition probabilities for the control group.
#'   - `Ps_treated`: Estimated transition probabilities for the treated group.
#'   - `Ps_control_empirical`: Estimated empirical transition probabilities for the control group.
#'   - `Ps_treated_empirical`: Estimated empirical transition probabilities for the treated group.
#'   - `pmfs_initial_control`: Estimated initial PMFs for the control group.
#'   - `pmfs_initial_treated`: Estimated initial PMFs for the treated group.
#'   - `pmfs_pre_treatment_control`: Estimated pre-treatment PMFs for the control group.
#'   - `pmfs_pre_treatment_treated`: Estimated pre-treatment PMFs for the treated group.
#'   - `treated_period`: The time period when treatment starts.
#'   - `Y`: The outcome names.
#'   The returned object is of class `"TransitionModel"`.
#'
#' @examples
#' set.seed(123)
#' scores <- list(
#'   priors = runif(3),
#'   priors_treated = runif(3),
#'   Ps_control = list(matrix(runif(9), ncol = 3), matrix(runif(9), ncol = 3)),
#'   Ps_treated = list(matrix(runif(9), ncol = 3), matrix(runif(9), ncol = 3)),
#'   Ps_control_empirical = list(matrix(runif(9), ncol = 3), matrix(runif(9), ncol = 3)),
#'   Ps_treated_empirical = list(matrix(runif(9), ncol = 3), matrix(runif(9), ncol = 3)),
#'   pmfs_initial_control = list(matrix(runif(9), ncol = 3)),
#'   pmfs_initial_treated = list(matrix(runif(9), ncol = 3)),
#'   pmfs_pre_treatment_control = list(matrix(runif(9), ncol = 3)),
#'   pmfs_pre_treatment_treated = list(matrix(runif(9), ncol = 3)),
#'   treated_period = 2,
#'   Y = list(names = c("Low", "Medium", "High"))
#' )
#' model <- scores_to_model(scores)
#' print(class(model))  # Should be "TransitionModel"
#'
#' @export
scores_to_model <- function(scores, weights = NULL) {
  # Extract data
  J <- length(scores$model$priors)
  Y <- scores$Y
  treated_period <- scores$treated_period
  T_max <- length(scores$Ps_control[[1]]) + 1
  prob_treated <- scores$prob_treated
  model <- scores$model

  # Split weights
  weights_treated <- NULL
  weights_control <- NULL
  if (!is.null(weights)) {
    if (length(weights) != nrow(scores$priors)) {
      stop("Length of weights must match number of observations.")
    }
    treated_indices <- scores$g != 0
    weights_treated <- weights[treated_indices]
    weights_control <- weights[!treated_indices]
  }

  # Compute priors
  priors <- score_to_pmf(scores$priors, model$priors, weights = weights)
  priors_treated <- score_to_pmf(scores$priors_treated, model$priors_treated, weights = weights_treated)

  # Compute transition matrices
  Ps_treated_empirical <- score_to_Ps(scores$Ps_treated_empirical, model$Ps_treated_empirical, weights = weights_treated)
  Ps_control_empirical <- score_to_Ps(scores$Ps_control_empirical, model$Ps_control_empirical, weights = weights_control)
  Ps_treated <- score_to_Ps(scores$Ps_treated, model$Ps_treated, weights = weights)
  Ps_control <- score_to_Ps(scores$Ps_control, model$Ps_control, weights = weights)

  # Note: Ps_control = Ps_treated is the same for treated and control in pre-treatment periods
  if (treated_period > 2) {
    Ps_treated <- lapply(1:J, function(j)
      c(Ps_treated[[j]][1:(treated_period-2)], Ps_treated_empirical[[j]][(treated_period - 1):(T_max - 1)]))
    Ps_control <- lapply(1:J, function(j)
      c(Ps_control[[j]][1:(treated_period-2)], Ps_control_empirical[[j]][(treated_period - 1):(T_max - 1)]))
  }

  # Compute PMFs
  pmfs_initial_control <- score_to_pmfs(scores$pmfs_initial_control, model$pmfs_initial_control,
                                        weights = weights_control)
  pmfs_initial_treated <- score_to_pmfs(scores$pmfs_initial_treated, model$pmfs_initial_treated,
                                        weights = weights_treated)
  pmfs_pre_treatment_control <- score_to_pmfs(scores$pmfs_pre_treatment_control, model$pmfs_pre_treatment_control,
                                              weights = weights_control)
  pmfs_pre_treatment_treated <- score_to_pmfs(scores$pmfs_pre_treatment_treated, model$pmfs_pre_treatment_treated,
                                              weights = weights_treated)

  # Return model
  model <- list(
    prob_treated = prob_treated,
    priors = priors,
    priors_treated = priors_treated,
    Ps_control = Ps_control,
    Ps_treated = Ps_treated,
    Ps_control_empirical = Ps_control_empirical,
    Ps_treated_empirical = Ps_treated_empirical,
    pmfs_initial_control = pmfs_initial_control,
    pmfs_initial_treated = pmfs_initial_treated,
    pmfs_pre_treatment_control = pmfs_pre_treatment_control,
    pmfs_pre_treatment_treated = pmfs_pre_treatment_treated,
    treated_period = treated_period,
    Y = Y
  )
  class(model) <- "TransitionModel"
  model
}
