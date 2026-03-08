#' @title Path Construction Functions for Empirical Examples
#'
#' @description
#' Functions for constructing file paths for plots, models, and other outputs
#' in empirical example analyses. These functions ensure consistent naming
#' conventions across the package.
#'
#' @name empirical-examples-paths
NULL

# Note: is_valid_varname() is defined in utility.R and exported from there

#' Get Time Variable Name
#'
#' Returns the name of the time variable used in the empirical specification.
#' Defaults to "time_index" if not explicitly provided in the specification.
#'
#' @param empirical_spec A list or object containing the empirical specification,
#'   including (optionally) the field \code{time_index_varname} representing
#'   the time variable.
#' @return A character string representing the name of the time variable.
#' @export
get_tname <- function(empirical_spec) {
  if (is_valid_varname(empirical_spec$time_index_varname)) {
    return(empirical_spec$time_index_varname)
  } else {
    return("time_index")
  }
}

#' Construct a Filename Suffix from Empirical Specification
#'
#' Generates a suffix string based on the clustering and cohort variables in
#' an empirical specification. This is useful for distinguishing output files
#' by modeling details such as clustered standard errors or cohort effects.
#'
#' @param empirical_spec A list or object containing the fields
#'   \code{cluster_varname} and \code{cohort_varname}, representing the names
#'   of the clustering and cohort variables (if any).
#' @param include_se_info A logical value indicating whether to include the
#'   standard error method in the suffix.
#' @return A character string beginning with a dash, followed by relevant
#'   specification details. Returns an empty string if no relevant fields
#'   are present.
#' @export
empirical_spec_suffix <- function(empirical_spec, include_se_info = TRUE) {
  cluster_varname <- empirical_spec$cluster_varname
  cohort_varname <- empirical_spec$cohort_varname
  sample_varname <- empirical_spec$sample_varname
  opts_se_method <- empirical_spec$opts_se_method

  suffix <- ""
  if (is_valid_varname(cluster_varname) && include_se_info) {
    suffix <- paste0(suffix, "-cluster_", cluster_varname)
  }
  if (is_valid_varname(cohort_varname)) {
    suffix <- paste0(suffix, "-cohort_", cohort_varname)
  }
  if (is_valid_varname(sample_varname)) {
    suffix <- paste0(suffix, "-sample_", sample_varname)
  }
  if (is_valid_varname(opts_se_method) && include_se_info) {
    suffix <- paste0(suffix, "-", opts_se_method)
  }

  suffix
}

# -----------------------------------------------------------------------------
# Directory Creation
# -----------------------------------------------------------------------------

#' Create Directory if Needed
#'
#' Ensures that a specified directory exists by creating it if it does not
#' already exist.
#'
#' @param directory A character string specifying the path to the directory.
#' @return None. Creates the directory if it does not exist.
#' @export
create_directory_if_needed <- function(directory) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
}

#' Construct Plot File Path
#'
#' Constructs a file path using the base directory, outcome names from the
#' empirical specification, and filename.
#'
#' @param base_dir A character string specifying the base directory.
#' @param empirical_spec A list containing the empirical specification details,
#'   including \code{outcome_varname} and \code{outcome_name}.
#' @param filename A character string for the filename.
#' @return A character string representing the full constructed path.
#' @export
construct_plot_file_path <- function(base_dir, empirical_spec, filename) {
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name
  paste0(base_dir, "/", outcome_varname, "/", outcome_name, "/", filename)
}

#' Build Plot Filename
#'
#' Generic helper function to build plot filenames with consistent naming
#' convention. This is used internally by the specific get_plot_*_directory
#' functions.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param plot_type A character string specifying the type of plot (e.g.,
#'   "averages", "ltatt", "transitions").
#' @param J An optional integer specifying the number of latent groups.
#'   Default is NULL.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @param include_se_info A logical value indicating whether to include SE info
#'   in spec suffix. Default is TRUE.
#' @param prefix_with_J A logical value indicating whether J should come before
#'   plot_type. Default is TRUE.
#' @return A character string representing the filename.
#' @keywords internal
build_plot_filename <- function(empirical_spec, plot_type, J = NULL,
                                suffix = "", include_se_info = TRUE,
                                prefix_with_J = TRUE) {
  outcome_varname <- empirical_spec$outcome_varname
  spec_suffix <- empirical_spec_suffix(empirical_spec,
                                       include_se_info = include_se_info)
  J_suffix <- if (!is.null(J)) paste0("-J", J) else ""

  if (prefix_with_J) {
    paste0("model-", outcome_varname, spec_suffix, J_suffix, "-", plot_type,
           suffix, ".pdf")
  } else {
    paste0("model-", outcome_varname, spec_suffix, "-", plot_type, J_suffix,
           suffix, ".pdf")
  }
}

#' Get Plot File Path
#'
#' Generic function to build plot file paths with consistent naming convention.
#' This can be used directly or through the specific get_plot_*_directory
#' wrapper functions.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param plot_type A character string specifying the type of plot (e.g.,
#'   "averages", "ltatt", "transitions").
#' @param J An optional integer specifying the number of latent groups.
#'   Default is NULL.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @param include_se_info A logical value indicating whether to include SE info
#'   in spec suffix. Default is TRUE.
#' @return A character string representing the full file path for the plot.
#' @export
get_plot_filepath <- function(empirical_spec, plot_type, J = NULL, suffix = "",
                              include_se_info = TRUE) {
  filename <- build_plot_filename(empirical_spec, plot_type, J, suffix,
                                  include_se_info)
  construct_plot_file_path("plots", empirical_spec, filename)
}

#' Construct Model File Path
#'
#' Constructs a file path for model files using the base directory, outcome
#' names, and filename.
#'
#' @param base_dir A character string specifying the base directory.
#' @param empirical_spec A list containing the empirical specification details.
#' @param filename A character string for the filename. Default is NULL.
#' @param subdir_with_outcome_name A logical value indicating whether to include
#'   the outcome name in the subdirectory. Default is FALSE.
#' @return A character string representing the full constructed path.
#' @export
construct_model_file_path <- function(base_dir, empirical_spec,
                                      filename = NULL,
                                      subdir_with_outcome_name = FALSE) {
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name

  if (is.null(outcome_varname)) {
    outcome_varname <- "outcome"
  }
  if (is.null(outcome_name) || !subdir_with_outcome_name) {
    return(paste0(base_dir, "/", outcome_varname, "/", filename))
  }
  paste0(base_dir, "/", outcome_varname, "/", outcome_name, "/", filename)
}

#' Generate Model Name
#'
#' Generates a standardized model name using the outcome variable name and
#' number of latent groups.
#'
#' @param J An integer specifying the number of latent groups.
#' @param empirical_spec A list containing the empirical specification details.
#' @return A character string representing the model file path.
#' @export
get_model_name <- function(J, empirical_spec) {
  outcome_varname <- empirical_spec$outcome_varname
  filename <- paste0("model-", outcome_varname,
                     empirical_spec_suffix(empirical_spec), "-J", J)
  construct_model_file_path("models", empirical_spec, filename = filename)
}

#' Generate DiD Model Name (TWFE)
#'
#' Generates a standardized model name for two-way fixed effects DiD models.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @return A character string representing the model file path.
#' @export
get_did_model_twfe_name <- function(empirical_spec) {
  outcome_varname <- empirical_spec$outcome_varname
  filename <- paste0("model-", outcome_varname,
                     empirical_spec_suffix(empirical_spec), "-did-twfe")
  construct_model_file_path("models", empirical_spec, filename = filename,
                            subdir_with_outcome_name = TRUE)
}

#' Generate DiD Model Name
#'
#' Generates a standardized model name for DiD models.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @return A character string representing the model file path.
#' @export
get_did_model_name <- function(empirical_spec) {
  outcome_varname <- empirical_spec$outcome_varname
  filename <- paste0("model-", outcome_varname,
                     empirical_spec_suffix(empirical_spec), "-did")
  construct_model_file_path("models", empirical_spec, filename = filename,
                            subdir_with_outcome_name = TRUE)
}

#' Generate Conditional Lag DiD Model Name
#'
#' Generates a standardized model name for conditional lag DiD models.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param lag An optional integer specifying the number of lags. Default is NULL.
#' @return A character string representing the model file path.
#' @export
get_did_conditional_lag_model_name <- function(empirical_spec, lag = NULL) {
  outcome_varname <- empirical_spec$outcome_varname
  filename <- paste0("model-", outcome_varname,
                     empirical_spec_suffix(empirical_spec),
                     "-did-conditional-lag-", lag)
  construct_model_file_path("models", empirical_spec, filename = filename,
                            subdir_with_outcome_name = TRUE)
}

# -----------------------------------------------------------------------------
# Plot Directory Functions
# -----------------------------------------------------------------------------

#' Get Plot Averages Directory
#'
#' Generates the file path for saving average plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param J An optional integer specifying the number of latent groups.
#'   Default is NULL.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for average plots.
#' @export
get_plot_averages_directory <- function(empirical_spec, J = NULL, suffix = "") {
  get_plot_filepath(empirical_spec, "averages", J = J, suffix = suffix,
                    include_se_info = FALSE)
}

#' Get Plot Averages (All Outcomes) Directory
#'
#' Generates the file path for saving average plots for all outcome states.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param J An optional integer specifying the number of latent groups.
#'   Default is NULL.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for average plots.
#' @export
get_plot_averages_all_outcome_directory <- function(empirical_spec, J = NULL,
                                                    suffix = "") {
  get_plot_filepath(empirical_spec, "averages-all-outcome", J = J,
                    suffix = suffix, include_se_info = FALSE)
}

#' Get Plot Conditional Averages Directory
#'
#' Generates the file path for saving conditional average plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for average plots.
#' @export
get_plot_conditional_averages_directory <- function(empirical_spec,
                                                    suffix = "") {
  get_plot_filepath(empirical_spec, "averages-conditional", suffix = suffix,
                    include_se_info = FALSE)
}

#' Get Plot LTATT Directory
#'
#' Generates the file path for saving LTATT (Latent Group Average Treatment
#' Effects) plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param J An integer specifying the number of latent groups.
#' @return A character string representing the file path for LTATT plots.
#' @export
get_plot_ltatt_directory <- function(empirical_spec, J) {
  get_plot_filepath(empirical_spec, "ltatt", J = J)
}

#' Get Plot LTATT by Pretreatment Outcomes Directory
#'
#' Generates the file path for saving LTATT plots stratified by pretreatment
#' outcomes.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param J An integer specifying the number of latent groups.
#' @return A character string representing the file path for LTATT plots.
#' @export
get_plot_ltatt_by_pretreatment_outcomes_directory <- function(empirical_spec,
                                                               J) {
  get_plot_filepath(empirical_spec, "ltatt-by-pretreatment-outcomes", J = J)
}

#' Get Plot Transitions Directory
#'
#' Generates the file path for saving transition probability plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param J An integer specifying the number of latent groups.
#' @param outcome_before An optional character string for the outcome state
#'   before transition. Default is NULL.
#' @param outcome_after An optional character string for the outcome state
#'   after transition. Default is NULL.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for transition plots.
#' @export
get_plot_transitions_directory <- function(empirical_spec, J,
                                           outcome_before = NULL,
                                           outcome_after = NULL,
                                           suffix = "") {
  transition_detail <- if (!is.null(outcome_before) &&
                           !is.null(outcome_after)) {
    paste0("-", outcome_before, "-to-", outcome_after)
  } else {
    ""
  }
  plot_type <- paste0("transitions", transition_detail)
  get_plot_filepath(empirical_spec, plot_type, J = J, suffix = suffix)
}

#' Get Plot Transitions Differences Directory
#'
#' Generates the file path for saving transition difference plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param J An integer specifying the number of latent groups.
#' @param outcome_before An optional character string for the outcome state
#'   before transition. Default is NULL.
#' @param outcome_after An optional character string for the outcome state
#'   after transition. Default is NULL.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for transition plots.
#' @export
get_plot_transitions_diffs_directory <- function(empirical_spec, J,
                                                  outcome_before = NULL,
                                                  outcome_after = NULL,
                                                  suffix = "") {
  transition_detail <- if (!is.null(outcome_before) &&
                           !is.null(outcome_after)) {
    paste0("-", outcome_before, "-to-", outcome_after)
  } else {
    ""
  }
  plot_type <- paste0("transitions-diffs", transition_detail)
  get_plot_filepath(empirical_spec, plot_type, J = J, suffix = suffix)
}

#' Get Plot ATT Directory
#'
#' Generates the file path for saving ATT (Average Treatment Effect on the
#' Treated) plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param J An integer specifying the number of latent groups.
#' @return A character string representing the file path for ATT plots.
#' @export
get_plot_att_directory <- function(empirical_spec, J) {
  get_plot_filepath(empirical_spec, "att", J = J)
}

#' Get Plot LTATT All Directory
#'
#' Generates the file path for saving combined LTATT plots for all cohorts.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @return A character string representing the file path for LTATT plots.
#' @export
get_plot_ltatt_all_directory <- function(empirical_spec) {
  get_plot_filepath(empirical_spec, "all-ltatt")
}

#' Get Plot ATT All Directory
#'
#' Generates the file path for saving combined ATT plots for all cohorts.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @return A character string representing the file path for ATT plots.
#' @export
get_plot_att_all_directory <- function(empirical_spec) {
  get_plot_filepath(empirical_spec, "all-att")
}

#' Get Event Study Plot Directory
#'
#' Generates the file path for saving event study plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for event study plots.
#' @export
get_plot_event_study_directory <- function(empirical_spec, suffix = "") {
  get_plot_filepath(empirical_spec, "event_study", suffix = suffix)
}

#' Get Plot ATT Comparison All Directory
#'
#' Generates the file path for saving ATT comparison plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @return A character string representing the file path for comparison plots.
#' @export
get_plot_att_comparison_all_directory <- function(empirical_spec) {
  get_plot_filepath(empirical_spec, "comparison-att")
}

#' Get Plot ATT from 0 Comparison All Directory
#'
#' Generates the file path for saving ATT from initial period comparison plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for comparison plots.
#' @export
get_plot_att_from_0_comparison_all_directory <- function(empirical_spec,
                                                          suffix = "") {
  get_plot_filepath(empirical_spec, "comparison-att_from_0", suffix = suffix)
}

#' Get Plot Counterfactuals Comparison Directory
#'
#' Generates the file path for saving counterfactual comparison plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for comparison plots.
#' @export
get_plot_counterfactuals_comparison_directory <- function(empirical_spec,
                                                           suffix = "") {
  get_plot_filepath(empirical_spec, "counterfactuals-comparison",
                    suffix = suffix)
}

#' Get Plot ATT Comparison Directory
#'
#' Generates the file path for saving ATT comparison plots.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for comparison plots.
#' @export
get_plot_att_comparison_directory <- function(empirical_spec, suffix = "") {
  get_plot_filepath(empirical_spec, "att-comparison", suffix = suffix)
}

#' Get Plot Counterfactuals by Pretreatment Outcomes Directory
#'
#' Generates the file path for saving counterfactual plots stratified by
#' pretreatment outcomes.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param J An integer specifying the number of latent groups. Default is NULL.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for counterfactual
#'   plots.
#' @export
get_plot_counterfactuals_by_pretreatment_outcomes_directory <- function(
    empirical_spec, J = NULL, suffix = "") {
  get_plot_filepath(empirical_spec, "counterfactuals-by-pretreatment-outcomes",
                    J = J, suffix = suffix)
}

#' Get Plot Counterfactuals by Outcomes Combined Directory
#'
#' Generates the file path for saving counterfactual plots stratified by both
#' initial and pretreatment outcomes.
#'
#' @param empirical_spec A list containing the empirical specification details.
#' @param J An integer specifying the number of latent groups. Default is NULL.
#' @param suffix An optional character string to add a suffix to the filename.
#'   Default is "".
#' @return A character string representing the file path for counterfactual
#'   plots.
#' @export
get_plot_counterfactuals_by_outcomes_combined_directory <- function(
    empirical_spec, J = NULL, suffix = "") {
  get_plot_filepath(empirical_spec, "counterfactuals-by-outcomes-combined",
                    J = J, suffix = suffix)
}
