
#' Estimate and Save DiD Model
#'
#' Estimates a did model using specified parameters and saves the result to a file.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param event_study_df A data frame containing the event study data.
#' @param opts A list of options for the model estimation.
#' @return The estimated model object.
#' @export
estimate_and_save_did_model <- function(empirical_spec, event_study_df, opts = NULL) {
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name
  outcomes_of_interest <- as.character(empirical_spec$outcomes_of_interest)
  outcomes_of_interest <- strsplit(outcomes_of_interest, ",")[[1]]
  cluster_varname <- NULL
  if (!is.na(empirical_spec$cluster_varname) &&
      empirical_spec$cluster_varname != "") {
    cluster_varname <- empirical_spec$cluster_varname
  }

  # Apply sample restriction
  event_study_df_restricted <- apply_sample_restriction(event_study_df, empirical_spec)

  yname <- outcome_varname
  tname <- get_tname(empirical_spec)
  idname <- "id"
  gname <- "treatment_cohort"

  # Perform DID
  did_df <- event_study_df_restricted %>%
    mutate("id_did" := as.numeric(factor(!!sym(idname)))) %>%
    mutate(outcome_value = !!sym(outcome_varname) %in% outcomes_of_interest) %>%
    mutate(outcome_value = outcome_value*1)
  did_model <- did::att_gt(yname = "outcome_value",
                           tname = tname,
                           idname = "id_did",
                           gname = gname,
                           clustervars = cluster_varname,
                           data = did_df)
  did_model$empirical_spec <- empirical_spec

  # Save the solved model
  model_name <- get_did_model_name(empirical_spec)
  save_model(did_model, model_name, paste0(model_name, ".RData"))

  did_model
}

#' Load or DiD Model
#'
#' Loads a saved model if it exists, otherwise estimates a new model and saves it.
#' If `force_estimate` is TRUE, a new model will be estimated and saved even if
#' an existing model file is found.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param event_study_df A data frame containing the event study data.
#' @param opts A list of options for the model estimation.
#' @param force_estimate A logical value. If TRUE, forces estimation even if a model file exists.
#' @return The loaded or estimated model object.
#' @export
load_or_estimate_did_model <- function(empirical_spec, event_study_df, opts = NULL,
                                       force_estimate = FALSE) {

  # Try loading the model
  model_name <- get_did_model_name(empirical_spec)
  filename <- paste0(model_name, ".RData")
  model <- load_model(filename)

  # Estimate and save the model if not loaded or force_estimate is TRUE
  if (force_estimate || is.null(model)) {
    model <- estimate_and_save_did_model(empirical_spec, event_study_df, opts)
  }

  # Re-estimate if saved model is incompatible with current did package version
  # (did >= 2.3.0 restructured DIDparams; patching individual fields is too brittle)
  if (!is.null(model) && is.null(model$DIDparams$time_periods_count)) {
    model <- estimate_and_save_did_model(empirical_spec, event_study_df, opts)
  }

  # Add in ATT aggregate
  model$att_aggregate <- mean(model$att[unique(model$group):length(model$att)])

  return(model)
}
