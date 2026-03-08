#' Load Model from File
#'
#' Loads a saved model from a specified file. If there are multiple objects in the file, an error is raised.
#'
#' @param filename A character string specifying the transition to the saved model file.
#' @return The loaded model object if only one object is present in the file. NULL if the model does not exist
#' @export
load_model <- function(filename) {
  # Return NULL if the file does not exist
  if (!file.exists(filename)) {
    return(NULL)
  }

  # Load the file and return the object if only one object is present
  e <- new.env()
  load(filename, envir = e)
  loaded <- as.list(e)

  if (length(loaded) == 1) {
    return(loaded[[1]])
  } else {
    stop(paste0("Multiple objects in model file: ", filename))
  }
  return(NULL)
}

#' Save Model to File
#'
#' This function saves an R object (typically a model) to a `.RData` file. The function assigns the object
#' to a variable in the global environment and saves it to the specified file name.
#'
#' @param model The R object to be saved, usually a trained model.
#' @param model_name A character string specifying the name to assign to the model in the global environment.
#' @param file_name A character string specifying the file name to save the model. Defaults to \code{model_name.RData}.
#' @param slim A logical; if TRUE (default), remove unnecessary components from the model
#'   (especially from bootstrap estimates) to reduce file size.
#' @param compress A character string specifying the compression type. Default is "xz" for
#'   best compression. Use "gzip" for faster but larger files, or FALSE for no compression.
#'
#' @return None. The function performs a side effect of saving the model to a file.
#' @export
save_model <- function(model, model_name, file_name = NULL, slim = TRUE, compress = "xz") {
  # Create the directory if needed
  create_directory_if_needed(dirname(file_name))

  # Save the model to the file
  if (is.null(file_name)) {
    file_name <- paste0(model_name, ".RData")
  }

  # Slim the model if requested
  if (slim) {
    model <- slim_model(model)
  }

  assign(model_name, model, envir = .GlobalEnv)
  save(list = model_name, file = file_name, compress = compress)
}

#' Slim Model for Storage
#'
#' Removes unnecessary components from a model to reduce storage size.
#' This is particularly useful for models with bootstrap estimates, where
#' each bootstrap sample stores large matrices (weights, posteriors) that
#' are not needed for computing confidence intervals.
#'
#' @param model A fitted transition model object.
#'
#' @return The model with unnecessary components removed.
#' @export
slim_model <- function(model) {

  # Components to remove from bootstrap estimates (large matrices not needed for CIs)
  bootstrap_remove <- c("weights", "posteriors", "unit_weights")

  # Components to remove from main model (can be recomputed if needed)
  # Note: Keep event_study_df as it's needed for plotting functions
  main_remove <- c("data")

  # Slim bootstrap estimates if present
  if (!is.null(model$bootstrap_estimates)) {
    model$bootstrap_estimates <- lapply(model$bootstrap_estimates, function(boot) {
      for (comp in bootstrap_remove) {
        boot[[comp]] <- NULL
      }
      boot
    })
  }

  # Remove large components from main model
  for (comp in main_remove) {
    model[[comp]] <- NULL
  }

  model
}

#' Estimate and Save DID Model
#'
#' Estimates a model using specified parameters and saves the result to a file.
#'
#' @param model_name A character string specifying the model name used for saving.
#' @param J An integer specifying the number of lags.
#' @param empirical_spec A list containing the empirical specification details.
#' @param event_study_df A data frame containing the event study data.
#' @param opts A list of options for the model estimation.
#' @return The estimated model object.
#' @export
estimate_and_save_model <- function(model_name, J,
                                    empirical_spec, event_study_df, opts = NULL) {
  # Extract parameters from empirical_spec
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name
  outcomes_of_interest <- empirical_spec$outcomes_of_interest
  cluster_varname <- NULL
  if (!is.na(empirical_spec$cluster_varname) &&
      empirical_spec$cluster_varname != "") {
    cluster_varname <- empirical_spec$cluster_varname
  }

  # Apply sample restriction
  event_study_df_restricted <- apply_sample_restriction(event_study_df, empirical_spec)

  # Apply options
  opts <- apply_opts_from_empirical_spec(opts, empirical_spec)

  # Set up the model parameters
  yname <- outcome_varname
  tname <- get_tname(empirical_spec)
  idname <- "id"
  gname <- "treatment_cohort"

  # Estimate the model
  solved <- ak::estimate_transition_model(
    yname = yname, gname = gname,
    tname = tname, idname = idname, cluster_varname = cluster_varname,
    outcomes_of_interest = outcomes_of_interest,
    event_study_df = event_study_df_restricted, J = J, opts = opts
  )

  # Add data as well for future reference (in case subsample is used)
  solved$event_study_df <- event_study_df_restricted
  solved$outcome_varname <- outcome_varname
  solved$outcome_name <- outcome_name

  # Save the solved model
  save_model(solved, model_name, paste0(model_name, ".RData"))

  solved
}

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
  tname <- "time_index"
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

#' Estimate and Save DiD Model
#'
#' Estimates a DiD model using Two-Way Fixed Effects (TWFE) with period-specific treatment effects and cluster robust standard errors.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param event_study_df A data frame containing the event study data.
#' @param opts A list of options for the model estimation.
#' @return The estimated model object with TWFE results and cluster robust standard errors.
#' @export
estimate_and_save_did_model_twfe <- function(empirical_spec, event_study_df, opts = NULL) {
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name
  outcomes_of_interest <- empirical_spec$outcomes_of_interest
  outcomes_of_interest <- strsplit(outcomes_of_interest, ",")[[1]]
  cluster_varname <- NULL
  if (!is.na(empirical_spec$cluster_varname) &&
      empirical_spec$cluster_varname != "") {
    cluster_varname <- empirical_spec$cluster_varname
  }

  # Apply sample restriction
  event_study_df_restricted <- apply_sample_restriction(event_study_df, empirical_spec)

  yname <- outcome_varname
  tname <- "time_index"
  idname <- "id"
  gname <- "treatment_cohort"

  # Prepare data for TWFE estimation
  did_df <- event_study_df_restricted %>%
    mutate("id_did" := as.numeric(factor(!!sym(idname)))) %>%
    mutate(outcome_value = !!sym(outcome_varname) %in% outcomes_of_interest) %>%
    mutate(outcome_value = outcome_value*1) %>%
    # Create treatment indicator for each post-treatment period
    mutate(treated = ifelse(!!sym(gname) <= !!sym(tname) & !!sym(gname) != 0, 1, 0)) %>%
    # Create period-specific treatment indicators
    mutate(periods_since_treatment = ifelse(treated == 1,
                                            !!sym(tname) - !!sym(gname),
                                            NA))

  # Get all unique post-treatment periods
  pre_treatment_periods <- 1:max(did_df$treatment_cohort)
  post_treatment_periods <- sort(unique(did_df$periods_since_treatment[!is.na(did_df$periods_since_treatment)]))

  # Create treatment dummies for each post-treatment period
  for (period in post_treatment_periods) {
    var_name <- paste0("treated_period_", period)
    did_df[[var_name]] <- ifelse(!is.na(did_df$periods_since_treatment) &
                                   did_df$periods_since_treatment == period, 1, 0)
  }

  # Build formula with all period-specific treatment effects
  treatment_vars <- paste0("treated_period_", post_treatment_periods)
  treatment_formula <- paste(treatment_vars, collapse = " + ")

  # Estimate TWFE model with period-specific effects using fixest
  if (!is.null(cluster_varname)) {
    # With clustering - fixest handles cluster robust SE automatically
    twfe_model <- fixest::feols(
      as.formula(paste0("outcome_value ~ ", treatment_formula, " | id_did + time_index")),
      data = did_df,
      cluster = cluster_varname
    )
  } else {
    # Without clustering - use heteroskedasticity robust SE
    twfe_model <- fixest::feols(
      as.formula(paste0("outcome_value ~ ", treatment_formula, " | id_did + time_index")),
      data = did_df,
      se = "hetero"
    )
  }

  # Extract results for all period-specific treatment effects
  coef_summary <- summary(twfe_model)

  # Extract coefficients safely
  all_coefs <- coef_summary$coeftable
  treatment_coefs <- all_coefs[rownames(all_coefs) %in% treatment_vars, , drop = FALSE]

  # Store results in did-style format for compatibility
  did_model <- list(
    # Period-specific treatment effects
    period_effects = list(
      periods = post_treatment_periods,
      att = treatment_coefs[, "Estimate"],
      se = treatment_coefs[, "Std. Error"],
      t = treatment_coefs[, "t value"],
      pval = treatment_coefs[, "Pr(>|t|)"],
      ci_lower = treatment_coefs[, "Estimate"] - 1.96 * treatment_coefs[, "Std. Error"],
      ci_upper = treatment_coefs[, "Estimate"] + 1.96 * treatment_coefs[, "Std. Error"]
    ),
    # Overall average treatment effect (mean of all period effects)
    att = (treatment_coefs[, "Estimate"]),
    se = treatment_coefs[, "Std. Error"],
    cluster_var = cluster_varname,
    n_obs = twfe_model$nobs,
    n_clusters = if (!is.null(cluster_varname)) length(unique(did_df[[cluster_varname]])) else NULL,
    method = "TWFE_fixest_period_specific"
  )

  # Add empirical specification
  did_model$empirical_spec <- empirical_spec

  # Add zeros in att in pretreatment periods like did model
  did_model$att <- c(rep(0, length(pre_treatment_periods)-1), did_model$att)
  did_model$se <- c(rep(0, length(pre_treatment_periods)-1), did_model$se)
  did_model$group <- max(did_df$treatment_cohort)

  # Save the solved model
  model_name <- get_did_model_twfe_name(empirical_spec)
  save_model(did_model, model_name, paste0(model_name, ".RData"))

  did_model
}

#' Load or Estimate Model
#'
#' Loads a saved model if it exists, otherwise estimates a new model and saves it.
#' If `force_estimate` is TRUE, a new model will be estimated and saved even if
#' an existing model file is found.
#'
#' @param J An integer specifying the number of lags.
#' @param empirical_spec A list containing the empirical specification details.
#' @param event_study_df A data frame containing the event study data.
#' @param opts A list of options for the model estimation.
#' @param force_estimate A logical value. If TRUE, forces estimation even if a model file exists.
#' @return The loaded or estimated model object.
#' @export
load_or_estimate_model <- function(J, empirical_spec, event_study_df, opts = NULL,
                                   force_estimate = FALSE) {
  # Try loading the model
  model_name <- get_model_name(J, empirical_spec)
  filename <- paste0(model_name, ".RData")
  model <- load_model(filename)

  # Estimate and save the model if not loaded or force_estimate is TRUE
  if (force_estimate || is.null(model)) {
    # Estimate and save the model if not loaded
    start_time <- proc.time()
    model <- estimate_and_save_model(model_name, J, empirical_spec,
                                     event_study_df, opts)
    elapsed <- (proc.time() - start_time)["elapsed"]
    cat(sprintf("Model J=%d estimated in %.1f seconds (%.1f minutes)\n", J, elapsed, elapsed / 60), file = stderr())
    gc()
  }

  # Update the outcome values to outcomes of interest;
  # Parse outcomes of interest by comma (,)
  outcomes_of_interest <- as.character(empirical_spec$outcomes_of_interest)
  outcomes_of_interest <- strsplit(outcomes_of_interest, ",")[[1]]
  model <- modify_outcomes_of_interest(model, outcomes_of_interest)

  # Update the model with the empirical specification
  model$empirical_spec <- empirical_spec

  return(model)
}

#' Load or DiD Model Using TWFE
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
load_or_estimate_did_model_twfe <- function(empirical_spec, event_study_df, opts = NULL,
                                   force_estimate = FALSE) {

  # Try loading the model
  model_name <- get_did_model_name(empirical_spec)
  filename <- paste0(model_name, ".RData")
  model <- load_model(filename)

  # Estimate and save the model if not loaded or force_estimate is TRUE
  if (force_estimate || is.null(model)) {
    # Estimate and save the model if not loaded
    return(estimate_and_save_did_model_twfe(empirical_spec, event_study_df, opts))
  }

  return(model)
}
