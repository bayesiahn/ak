#' Perform Empirical Analysis
#'
#' This function performs an empirical analysis based on a user-defined empirical specification.
#' It estimates a set of transition models and a difference-in-differences (DiD) model
#' using the supplied event study data.
#'
#' @param empirical_spec A list containing empirical specifications, including `J_max` which
#'   determines the number of transition models to estimate.
#' @param event_study_df A data frame containing the data required for the event study analysis.
#' @param opts Optional list of additional arguments or configuration options used by the model
#'   estimation functions.
#' @param force_estimate A logical value. If TRUE, forces re-estimation even if saved model files exist.
#'
#' @return
#' A list or data frame containing the results of the empirical analysis, such as
#' model summaries, selected paths, or evaluation metrics.
#'
#' @export
perform_empirical_analysis <- function(empirical_spec, event_study_df, opts = NULL,
                                       force_estimate = FALSE) {
  J_max <- empirical_spec$J_max
  J1_lag_max <- empirical_spec$J1_lag_max

  # Load or estimate transition models
  transition_models <- lapply(1:J_max, load_or_estimate_model,
                              empirical_spec = empirical_spec,
                              event_study_df = event_study_df,
                              opts = opts,
                              force_estimate = force_estimate)

  # Load or estimate DiD model
  did_result <- load_or_estimate_did_model(empirical_spec, event_study_df, opts,
                                           force_estimate = force_estimate)

  # Load or estimate more lag model
  did_conditional_lag_models <- lapply(1:J1_lag_max, function(lag)
                            load_or_estimate_did_conditional_lag_model(
                              empirical_spec = empirical_spec,
                              event_study_df = event_study_df,
                              opts = opts, lag = lag,
                              force_estimate = force_estimate))
  did_conditional_lag_model <- did_conditional_lag_models[[1]]

  # Return them as a list
  empirical_analysis_result <- list(empirical_spec = empirical_spec,
       transition_models = transition_models,
       did_result = did_result,
       did_conditional_lag_models = did_conditional_lag_models,
       did_conditional_lag_model = did_conditional_lag_model)

  class(empirical_analysis_result) <- c("empirical_analysis_result", class(empirical_analysis_result))
  empirical_analysis_result
}
