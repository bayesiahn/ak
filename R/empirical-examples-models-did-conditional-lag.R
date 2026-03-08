
#' Estimate and Save DiD Model
#'
#' Estimates a did model using specified parameters and saves the result to a file.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param event_study_df A data frame containing the event study data.
#' @param lag An integer specifying the number of lags to include in the model.
#' @param opts A list of options for the model estimation.
#' @return The estimated model object.
#' @export
estimate_and_save_did_conditional_lag_model <- function(empirical_spec, event_study_df, lag = 2, opts = NULL) {
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
  did_df_all <- event_study_df_restricted %>%
    mutate(id := as.numeric(factor(!!sym(idname)))) %>%
    mutate(y = !!sym(outcome_varname),
          y_numeric = 1*(!!sym(outcome_varname) %in% outcomes_of_interest),
          t = !!sym(tname)) %>%
    attach_pretreatment_outcome_sequence(lag = lag) %>%
    attach_pretreatment_outcome_sequence(lag = 1, "past_outcome")

  trim_results <- trim_rare_sequences(did_df_all, 5)
  did_df <- trim_results$did_df_trimmed
  did_df$past_outcome_sequence <- droplevels(did_df$past_outcome_sequence)

  did_model <- did::att_gt(yname = "y_numeric",
                           tname = "t",
                           idname = "id",
                           gname = gname,
                           clustervars = cluster_varname,
                           data = did_df,
                           xformla = ~ past_outcome_sequence,
                           est_method = "reg")
  did_model$empirical_spec <- empirical_spec
  did_model$trim_summary <- trim_results$summary
  did_model$data_untrimmed <- did_df_all

  # Save the solved model
  model_name <- get_did_conditional_lag_model_name(empirical_spec, lag = lag)
  save_model(did_model, model_name, paste0(model_name, ".RData"))

  did_model
}

#' Attach Encoded Sequence of Pretreatment Outcome(s) to Panel Data
#'
#' Constructs a past outcome sequence using a one-hot encoding of the outcome variable in periods
#' prior to treatment, and attaches the encoded sequence as a categorical covariate to the input panel data.
#'
#' @param did_df A data frame in long format with at least the following columns:
#'   \itemize{
#'     \item{\code{id}}: Unit identifier.
#'     \item{\code{t}}: Time period index.
#'     \item{\code{y}}: Outcome variable (can be binary or categorical).
#'     \item{\code{treatment_cohort}}: Group-level treatment assignment time.
#'   }
#' @param lag An integer specifying the number of past time periods to consider prior to treatment.
#'   Defaults to 2.
#' @param varname A string specifying the name of the new column to be added to \code{did_df}.
#'
#' @return A modified version of \code{did_df} with an additional column:
#' \describe{
#'   \item{past_outcome_sequence}{A factor encoding of the concatenated one-hot past outcome sequence
#'   (e.g., "y_t-2_0_y_t-1_1").}
#' }
#'
#' @details
#' The function filters for time periods strictly before the treatment period and within a lag window.
#' Each unique outcome-time combination is one-hot encoded, and the result is collapsed into a string to
#' serve as a covariate capturing pre-treatment history. This covariate can be used in DiD regressions
#' as a flexible way to adjust for past outcomes.
#'
#' @seealso [load_or_estimate_did_conditional_lag_model()], [did::att_gt()]
attach_pretreatment_outcome_sequence <- function(did_df, lag = 2, varname = "past_outcome_sequence") {
  # Extract treated period
  treated_period <- max(did_df$treatment_cohort)

  # Filter out unrelated far-past outcomes and post-treatment outcomes
  outcome_sequence_df <- did_df %>%
    filter(t < treated_period) %>%
    filter(t > treated_period - lag - 1) %>%
    select(all_of(c("id", "t", "y")))

  # Create outcome-time column
  outcome_sequence_df <- outcome_sequence_df %>%
    mutate(y = (as.character(y))) %>%
    mutate(outcome_time = paste0("y_t", t, "_", y))

  # Encode the outcome sequence
  encoded <- outcome_sequence_df %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(
      id_cols = c("id"),
      names_from = outcome_time,
      values_from = y,
      values_fill = ""
    ) %>%
    select(order(names(.)))

  # Generate past outcome sequence
  past_outcome_sequence <- do.call(paste, c(encoded %>% select(-id), sep = "-"))

  # Derive the actual sequence by collapsing non-empty indicators
  outcome_sequence_df <- encoded %>%
    tidyr::pivot_longer(
      cols = -id,
      names_to = "outcome_time",
      values_to = "value"
    ) %>%
    filter(value != "") %>%
    tidyr::separate(outcome_time, into = c("prefix", "t", "status"), sep = "_", remove = FALSE) %>%
    select(id, t, status) %>%
    tidyr::pivot_wider(
      names_from = t,
      values_from = status,
      names_prefix = "t",
      names_sort = TRUE
    ) %>%
    mutate(past_outcome_sequence_added = do.call(paste, c(across(starts_with("t")), sep = "-"))) %>%
    select(id, past_outcome_sequence_added)

  # Return attached df
  did_df %>%
    left_join(outcome_sequence_df, by = "id") %>%
    mutate(!!sym(varname) := as.factor(past_outcome_sequence_added)) %>%
    select(-past_outcome_sequence_added)
}

#' Load or estimate a model with more lags, using did package
#'
#' Loads a saved model if it exists, otherwise estimates a new model and saves it.
#' If `force_estimate` is TRUE, a new model will be estimated and saved even if
#' an existing model file is found.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param event_study_df A data frame containing the event study data.
#' @param lag An integer specifying the number of lags to include in the model.
#' @param opts A list of options for the model estimation.
#' @param force_estimate A logical value. If TRUE, forces estimation even if a model file exists.
#' @return The loaded or estimated model object.
#' @export
load_or_estimate_did_conditional_lag_model <- function(empirical_spec, event_study_df, lag = 2, opts = NULL,
                                       force_estimate = FALSE) {
  # Try loading the model
  model_name <- get_did_conditional_lag_model_name(empirical_spec, lag = lag)
  filename <- paste0(model_name, ".RData")
  model <- load_model(filename)

  # Estimate and save the model if not loaded or force_estimate is TRUE
  if (force_estimate || is.null(model)) {
    # Estimate and save the model if not loaded
    model <- estimate_and_save_did_conditional_lag_model(empirical_spec, event_study_df, lag, opts)
  }

  # Patch for compatibility with did package >= 2.3.0: add faster_mode if missing
  if (is.null(model$DIDparams$faster_mode)) {
    model$DIDparams$faster_mode <- TRUE
  }

  # Patch for compatibility with did package >= 2.3.0: convert data to data.table
  if (!inherits(model$DIDparams$data, "data.table")) {
    model$DIDparams$data <- data.table::as.data.table(model$DIDparams$data)
  }

  # Add in ATT aggregate and other info
  model$att_aggregate <- mean(model$att[unique(model$group):length(model$att)])
  model$lag <- lag

  return(model)
}

#' Trim rare past outcome sequences by control group prevalence
#'
#' This function trims rows from a DID-style dataframe where a `past_outcome_sequence` is
#' uncommon in the control group (i.e., occurs fewer than `n_trim_threshold` times at baseline
#' period, defined as `t == 1`), but is present in the treated group. This helps ensure that
#' each treated sequence has sufficient support in the control group for valid comparisons.
#'
#' @param did_df A data frame that includes at least the columns `treatment_cohort`,
#'   `past_outcome_sequence`, and `t`.
#' @param n_trim_threshold An integer. A sequence will be trimmed if it occurs fewer than
#'   this number of times in the control group but is present in the treated group.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{did_df_trimmed}{The data frame after removing rows with rare sequences.}
#'   \item{summary}{A list with trimming metadata, including:}
#'   \describe{
#'     \item{past_outcome_sequence_table}{A wide-format summary table showing sequence counts
#'     for treated and control groups at baseline.}
#'     \item{trimmed_sequences}{A character vector of sequences that were trimmed.}
#'     \item{n_trim_threshold}{The threshold used to define rarity.}
#'   }
#' }
#'
#' @export
trim_rare_sequences <- function(did_df, n_trim_threshold = 5) {
  past_outcome_sequence_table <- did_df %>%
    filter(t == 1) %>%
    mutate(treated = ifelse(treatment_cohort > 0, "treated", "control")) %>%
    group_by(treated, past_outcome_sequence) %>%
    summarize(n = n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = treated,
      values_from = n,
      names_prefix = "n_",
      values_fill = 0
    )

  # Identify rare sequences
  rare_seqs <- past_outcome_sequence_table %>%
    filter(n_control < n_trim_threshold & n_treated > 0) %>%
    select(past_outcome_sequence)

  if (length(rare_seqs$past_outcome_sequence) > 0) {
    message("Trimming rare past outcome sequences: ",
            paste(rare_seqs$past_outcome_sequence, collapse = ", "),
            " (occurs fewer than ", n_trim_threshold,
            " times in control units, while observed among at least one treated units.)")
    message("Table:",
            "\n", paste(capture.output(print(past_outcome_sequence_table)), collapse = "\n"),
            "\nRare sequences: ", paste(rare_seqs$past_outcome_sequence, collapse = ", "),
            "\nTrim threshold: ", n_trim_threshold, "\n",
            "\nTable:\n")
    print(past_outcome_sequence_table)
  }

  did_df_trimmed <- did_df %>%
    filter(!past_outcome_sequence %in% rare_seqs$past_outcome_sequence)

  trim_info <- list(
    past_outcome_sequence_table = past_outcome_sequence_table,
    trimmed_sequences = rare_seqs$past_outcome_sequence,
    n_trim_threshold = n_trim_threshold
  )

  list(did_df_trimmed = did_df_trimmed,
       summary = trim_info)
}
