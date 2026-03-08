#' Estimate a Transition Model for Treatment Effects Analysis
#'
#' This function estimates a transition model for treatment effects analysis by
#' executing both the first and second stages of the algorithm. It first
#' constructs a transition model based on pretreatment data and then extends it to
#' include post-treatment periods.
#'
#' @param J The number of latent groups in the model.
#' @param y A matrix or data frame of outcomes for each unit over time.
#' @param g A vector indicating the treatment status for each unit.
#' @param c A vector indicating the cluster for each unit (optional)
#' @param outcomes_of_interest A string or numeric value for the outcome of interest. Default is maximum
#' @param unit_weights A numeric vector of weights for each unit, or \code{NULL}.
#'        If \code{NULL}, equal weights are assumed.
#' @param transition_model_guesses A list of initial transition model guesses.
#'        If \code{NULL}, initial guesses are generated internally.
#' @param opts A list of options for the algorithm.
#'
#' @return A list containing the estimated transition model attached with
#' weights, unit_weights, posteriors, with Ps_control_empirical and Ps_treated_empirical
#'
#' @export
estimate_transition_model_from_wide_matrix <- function(J, y, g, c = NULL,
                                                       outcomes_of_interest = NULL,
                                                       unit_weights = NULL,
                                                 transition_model_guesses = NULL,
                                                 opts = NULL) {
  # Filter out rows with any NA values in y
  # Panel data may have missing observations which become NA in wide format
  complete_rows <- complete.cases(y)
  if (sum(!complete_rows) > 0) {
    message(sprintf("Filtering out %d rows with missing values (keeping %d complete cases)",
                    sum(!complete_rows), sum(complete_rows)))
    y <- y[complete_rows, , drop = FALSE]
    g <- g[complete_rows]
    if (!is.null(c)) {
      c <- c[complete_rows]
    }
  }

  # Validate and handle unit weights
  unit_weights <- check_unit_weights_validity(y, unit_weights)

  # Prepare data - exclude NA from unique values
  Y_values <- sort(unique(c(y[!is.na(y)])))
  treated_period <- max(g, na.rm = TRUE)
  y_indices_matrix <- y_to_y_indices(y, Y_values)
  data_for_est <- split_data_for_estimation(y_indices_matrix, g)

  # First stage
  first_stage_result <- first_stage(data_for_est, g, J,
                                    unit_weights = unit_weights,
                                    transition_model_guesses = transition_model_guesses,
                                    opts = opts)
  transition_model <- first_stage_result$transition_model
  transition_model$Y <- generate_Y_t(Y_values, outcomes_of_interest)
  posteriors <- first_stage_result$posteriors

  # Attach post-EM parameters
  post_em_params <- estimate_post_em_params(posteriors, data_for_est, unit_weights, opts)
  transition_model_estimate <- c(transition_model, post_em_params)

  if (!check_if_simplify(opts)) {
    transition_model_estimate <- c(
      transition_model_estimate,
      list(data = list(y = y, g = g, c = c))
    )
  }

  class(transition_model_estimate) <- "TransitionModel"

  if (get_se_method(opts) == "bootstrap") {
    # TODO: merge multiplier bootstrap accordingly
    opts_bootstrap <- opts
    opts_bootstrap$se_method <- "none"
    opts_bootstrap$simplify <- TRUE
    bootstrap_estimates <- bootstrap_transition_models_from_wide_matrix(J, y, g, c,
                                                                        outcomes_of_interest = outcomes_of_interest,
      transition_model_guesses = list(transition_model), opts = opts_bootstrap)
    transition_model_estimate$bootstrap_estimates <- bootstrap_estimates
  }

  if (get_se_method(opts) == "bootstrap_multiplier") {
    if (J > 1) {
      stop("Bootstrap multiplier method is not implemented for J > 1")
    }
    bootstrap_results <- multiplier_bootstrap(y, g, transition_model_estimate, opts = opts)
    transition_model_estimate$bootstrap_estimates <- bootstrap_results$bootstrap_estimates
    transition_model_estimate$scores <- bootstrap_results$scores
    transition_model_estimate$bootstrap_multipliers <- bootstrap_results$bootstrap_multipliers
  }

  return(reorder_transition_model(transition_model_estimate))
}

#' Estimate a Transition Model for Treatment Effects Analysis
#'
#' This function estimates a transition model for treatment effects analysis by
#' executing both the first and second stages of the algorithm. It first
#' constructs a transition model based on pretreatment data and then extends it to
#' include post-treatment periods.
#'
#' @param J The number of latent transition groups in the model.
#' @param event_study_df long panel data
#' @param yname name of dependent variable
#' @param gname name of group membership variable
#' @param idname name of id variable
#' @param tname name of time variable
#' @param cluster_varname name of cluster variable (optional)
#' @param opts A list of options for the algorithm.
#'
#' @return A list containing the estimated transition model, weights, posteriors, weights, and LTATTs.
#'
#' @export
estimate_transition_model <- function (J, event_study_df, yname, gname, idname, tname, cluster_varname = NULL,
                                       outcomes_of_interest = NULL,
                                       opts = list()) {
  # get y and gs
  transformed <- long_df_to_y_and_g(event_study_df, yname, gname, idname, tname, cluster_varname)

  y <- transformed$y
  g <- transformed$g
  c <- transformed$c
  return(estimate_transition_model_from_wide_matrix(J, y, g, c,
                                                    outcomes_of_interest = outcomes_of_interest,
                                                    opts = opts))
}

#' Check if hard classification is performed
#'
#' Checks if hard classification is performed in the first stage of the transition model.
#'
#' @param opts A list of options for the transition model.
#' @return A boolean indicating whether hard classification is performed; default is FALSE.
check_if_hard_classification <- function(opts = list()) {
  if (is.null(opts) || !is.list(opts)) {
    return(FALSE)
  }
  if (is.null(opts$hard_classification)) {
    return(FALSE)
  }
  return(opts$hard_classification)
}

#' Perform hard classification
#'
#' Performs hard classification on the posteriors from the first stage of the transition model.
#'
#' @param posteriors A matrix of posteriors from the first stage of the transition model.
#' @return A matrix of weights for the second stage of the transition model.
#' @export
hard_classification <- function(posteriors) {
  # Generate weights matrix
  weights <- matrix(0, ncol = ncol(posteriors), nrow = nrow(posteriors))

  # Find the indices of the max values in each row
  max_indices <- max.col(posteriors, ties.method = "first")

  # Convert these indices into a linear index
  linear_indices <- cbind(seq_len(nrow(posteriors)), max_indices)

  # Set these indices in the weights matrix to TRUE (or 1)
  weights[linear_indices] <- 1

  return(weights)
}

#' Check if returned estimated model contains minimal elements
#'
#' Checks if returned estimated model contains minimal elements after estimate_model
#'
#' @param opts A list of options for the transition model.
#' @return A boolean indicating whether simplified estimate is returned; default is FALSE.
check_if_simplify <- function(opts = list()) {
  if (is.null(opts) || !is.list(opts)) {
    return(FALSE)
  }
  if (is.null(opts$simplify)) {
    return(FALSE)
  }
  return(opts$simplify)
}
