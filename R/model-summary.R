#' Summary Method for TransitionModel Objects
#'
#' This function generates a summary for objects of class \code{TransitionModel}.
#'
#' @param model An object of class \code{TransitionModel}.
#' @param ... Additional arguments passed to or from other methods (not used here).
#' @return A list containing the summary elements of the transition model.
#' @export
summarize_transition_model <- function(model, ...) {
  # Extract relevant model components
  Y_values <- model$Y$values
  J <- length(model$Ps_control)
  treated_period <- model$treated_period
  T_max <- length(model$Ps_control[[1]]) + 1

  # Compute ATTs; if Ps_control_from_pre and Ps_treated_from_pre are NULL,
  # compute LTATTs under the assumption of Markovian transitions.
  LTATT <- list()
  if (is.null(model$Ps_control_from_pre) || is.null(model$Ps_treated_from_pre)) {
    LTATTs <- lapply(1:J, function(j) {
      compute_LTATTs_under_markovian(model$pmfs_pre_treatment_treated[[j]],
                                     model$Ps_treated[[j]][(treated_period-1):(T_max-1)],
                                     model$Ps_control[[j]][(treated_period-1):(T_max-1)],
                                     Y_values)
    })
  } else {
    LTATTs <- lapply(1:J, function(j) {
      compute_LTATTs(model$pmfs_pre_treatment_treated[[j]],
                     model$Ps_treated_from_pre[[j]][(treated_period-1):(T_max-1)],
                     model$Ps_control_from_pre[[j]][(treated_period-1):(T_max-1)],
                     Y_values)
    })
  }

  ATT <- aggregate_ATT_by_latent_type(LTATTs, model$priors_treated)

  # Compute expected outcomes
  expected_outcomes_empirical_df <- get_expected_outcomes_df(model)

  # Compute transition probability differences
  Ps_diffs <- lapply(1:J, function(j) {
    Map(`-`,model$Ps_treated[[j]], model$Ps_control[[j]])
  })
  Ps_empirical_diffs <- lapply(1:J, function(j) {
    Map(`-`,model$Ps_treated_empirical[[j]], model$Ps_control_empirical[[j]])
  })

  # Compute AIC and BIC if data is available
  aic <- NA
  bic <- NA
  if (!is.null(model$data)) {
    n <- length(model$data$g)
    ll <- model$log_likelihood
    k <- length(vectorize_transition_model(model))
    aic <- 2 * k - 2 * ll
    bic <- log(n) * k - 2 * ll
  }

  return(list(LTATTs = LTATTs,
              ATT = ATT,
              Ps_diffs = Ps_diffs,
              Ps_empirical_diffs = Ps_empirical_diffs,
              aic = aic,
              bic = bic,
              model = model))
}

#' Summarize Multiple Transition Models
#'
#' This function generates summaries for a list of \code{TransitionModel} objects.
#'
#' @param models A list of \code{TransitionModel} objects.
#' @return A list of summaries for each model in \code{models}, where each summary is a \code{TransitionModel} object.
#' @export
summarize_transition_models <- function(models) {
  # Use mclapply to summarize each model in parallel
  summaries <- parallel::mclapply(models, summarize_transition_model)

  return(summaries)
}



#' @method print TransitionModel
#' @export
print.TransitionModel <- function(x) {
  # Extract parameters
  priors <- x$priors
  priors_treated <- x$priors_treated
  treated_period <- x$treated_period
  model_summary <- summarize_transition_model(x)
  ATT <- model_summary$ATT
  LTATTs_matrix <- do.call(cbind,model_summary$LTATTs)

  # Extract names
  J <- length(priors)
  post_treatment_periods <- treated_period:(treated_period+length(ATT)-1)

  # Name estimates
  names(priors) <- paste("j =", 1:J)
  names(priors_treated) <- paste("j =", 1:J)
  names(ATT) <- paste("t =", post_treatment_periods)
  colnames(LTATTs_matrix) <- paste("j =", 1:J)
  rownames(LTATTs_matrix) <- paste("t =", post_treatment_periods)

  # Print info
  if (!is.null(x$empirical_spec)) {
    cat("Model for", empirical_spec_to_plot_title(x$empirical_spec), "\n")
  }
  cat("Transition Model Estimate:\n")
  cat("Number of Latent Groups (J):", J, "\n")
  cat("Prior probabilities P(Z_i = j): \n")
  print(round(priors,3))
  cat("Prior probabilities for the treated P(Z_i = j | D_i = 1): \n")
  print(round(priors_treated,3))
  cat("Average Treatment Effects (ATT):\n")
  print(round(mean(ATT),3))
  cat("Average Treatment Effects (ATT) by time:\n")
  print(round(ATT,3))
  cat("Latent Type Average Treatment Effects (LTATTs) by time and group:\n")
  print(round(LTATTs_matrix,3))

  # Print AIC and BIC if available
  if (is_valid_value(model_summary$aic)) {
    cat("AIC:", round(model_summary$aic, 2), "\n")
  }
  if (is_valid_value(model_summary$bic)) {
    cat("BIC:", round(model_summary$bic, 2), "\n")
  }
}

