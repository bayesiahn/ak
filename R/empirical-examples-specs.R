
#' Create Plot Title from Empirical Specification
#'
#' Constructs a plot title using the outcome name from an empirical specification,
#' optionally appending an annotation in parentheses.
#'
#' @param empirical_spec A list or object containing at least the field `outcome_name`,
#'   typically representing the outcome variable of interest.
#' @param prefix_title Optional string to be appended in before the title.
#' @param display_se_info Logical. If TRUE, includes cohort information in the title.
#'
#' @return A character string suitable for use as a plot title.
#'
#' @examples
#' spec <- list(outcome_name = "Employment Rate")
#' empirical_spec_to_plot_title(spec)
#' empirical_spec_to_plot_title(spec, "Placebo test")
#' @export
empirical_spec_to_plot_title <- function(empirical_spec, prefix_title = NULL,
                                         display_se_info = TRUE) {
  # Use outcome name and suffix to create a title for the plot
  title <- empirical_spec$outcome_name
  title <- ifelse(is.null(title), "", title)
  cluster_varname <- empirical_spec$cluster_varname
  cluster_name <- empirical_spec$cluster_name
  cohort_varname <- empirical_spec$cohort_varname
  cohort_name <- empirical_spec$cohort_name
  sample_varname <- empirical_spec$sample_varname
  sample_name <- empirical_spec$sample_name
  opts_se_method <- empirical_spec$opts_se_method
  opts_se_method_name <- empirical_spec$opts_se_method_name

  # Add prefix_title if needed
  if (!is.null(prefix_title)) {
    title <- ifelse(is.null(title) || title == "", title,
                    paste0(prefix_title, ": ", title))
  }
  if (is_valid_varname(cluster_varname) && display_se_info) {
    title <- paste0(title, ", Clustered by ",
                    ifelse(is.null(cluster_name) || cluster_name == "",
                           cluster_varname, cluster_name))
  }
  if (is_valid_varname(cohort_varname)) {
    title <- paste0(title, " (",
                    ifelse(is.null(cohort_name) || cohort_name == "",
                           cohort_varname, cohort_name), ")")
  }
  if (is_valid_varname(sample_varname) && is_valid_varname(sample_name)) {
    title <- paste0(title, ", Sample: ",
                    ifelse(is.null(sample_name) || sample_name == "",
                           sample_varname, sample_name))
  }
  if (is_valid_varname(opts_se_method) && display_se_info) {
    title <- paste0(title, ", ",
                    ifelse(is.null(opts_se_method_name) || opts_se_method_name == "",
                           opts_se_method, opts_se_method_name))
  }

  title
}

#' Print method for empirical_spec objects
#'
#' Prints a human-readable summary of the empirical specification object.
#'
#' @param x An object of class \code{"empirical_spec"}.
#' @param ... Additional arguments (ignored).
#' @export
print.empirical_spec <- function(x, ...) {
  cat("Empirical Specification: ", empirical_spec_to_plot_title(x), "\n", sep = "")
  cat("Suffix for plots: ", empirical_spec_suffix(x), "\n", sep = "")

  # Optionally print core fields
  if (is_valid_varname(x$outcome_varname)) {
    cat("Outcome Variable: ", x$outcome_varname, "\n", sep = "")
  }
  if (!is.null(x$outcomes_of_interest)) {
    cat("Outcomes of Interest: ", x$outcomes_of_interest, "\n", sep = "")
  }
  if (!is.null(x$J_max)) {
    cat("Max J (number of latent groups): ", x$J_max, "\n", sep = "")
  }
  if (!is.null(x$J1_lag_max)) {
    cat("Max lag (lags for J = 1 models): ", x$J1_lag_max, "\n", sep = "")
  }
  if (is_valid_varname(x$cluster_varname)) {
    cat("Clustered by: ", x$cluster_varname, "\n", sep = "")
  }
  if (is_valid_varname(x$cohort_varname)) {
    cat("Status Variable: ", x$cohort_varname, "\n", sep = "")
  }
  if (is_valid_varname(x$sample_varname)) {
    cat("Subsample of: ", x$sample_varname, "\n", sep = "")
  }
  if (is_valid_varname(x$opts_se_method)) {
    cat("Standard Error Method: ", x$opts_se_method, "\n", sep = "")
  } else {
    cat("(Standard Error Method: Not specified; default is weighted bootstrap)\n")
  }

  if (!is_valid_varname(x$cluster_varname)) {
    cat("(Cluster: No clustering)\n")
  }
  if (!is_valid_varname(x$sample_varname)) {
    cat("(Sample: Entire population)\n")
  }

  invisible(x)
}

#' Convert a Data Frame of Empirical Specifications to a List of empirical_spec Objects
#'
#' This function takes a data frame of empirical specifications and converts each row into a list,
#' assigning the class \code{"empirical_spec"} to each resulting list. The output is a list of these
#' empirical specification objects.
#'
#' @param empirical_specs_df A data frame where each row represents an empirical specification.
#'
#' @return A list of objects of class \code{"empirical_spec"}, each corresponding to a row in the input data frame.
#'
#' @importFrom purrr pmap map
#' @export
empirical_specs_df_to_empirical_specs <- function(empirical_specs_df) {
  empirical_specs_df %>%
    purrr::pmap(~ list(...)) %>%
    purrr::map(~ `class<-`(.x, c("empirical_spec", class(.x))))
}

#' Apply Sample Restrictions to Event Study Data
#'
#' Prepares an event study data frame by selecting relevant variables, handling
#' missing outcomes, filtering by a sample indicator, and renaming the cohort
#' variable if specified in the empirical specification.
#'
#' @param event_study_df A data frame containing event study variables,
#'   including unit identifiers, time indices, treatment cohorts, outcome
#'   values, and optional cluster or sample indicators.
#' @param empirical_spec A list or object containing at least the fields:
#'   \itemize{
#'     \item \code{outcome_varname}: Name of the outcome variable.
#'     \item \code{outcome_name}: Descriptive label for the outcome.
#'     \item \code{cluster_varname}: Optional name of the clustering variable.
#'     \item \code{sample_varname}: Optional binary indicator for sample
#'       inclusion.
#'     \item \code{cohort_varname}: Optional name of the treatment cohort
#'       variable.
#'   }
#'
#' @return A modified data frame filtered to include only relevant observations,
#'   with missing outcomes set to 2 and cohort variable renamed to
#'   \code{treatment_cohort} if applicable.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Selects only the columns relevant for analysis based on the
#'     empirical specification.
#'   \item Replaces \code{NA} values in the outcome variable with \code{2}.
#'   \item Filters the data to only include rows where \code{sample_varname == 1},
#'     if applicable.
#'   \item Renames the cohort variable to \code{treatment_cohort}, if specified.
#' }
#'
#' @export
apply_sample_restriction <- function(event_study_df, empirical_spec) {
  # Extract relevant variables from empirical_spec
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name
  cluster_varname <- empirical_spec$cluster_varname
  cohort_varname <- empirical_spec$cohort_varname
  sample_varname <- empirical_spec$sample_varname
  period_varname <- get_tname(empirical_spec)
  covariate_table_varnames <- empirical_spec$covariate_table_varnames
  if (is.na(covariate_table_varnames)) {
    covariate_table_varnames <- ""
  }
  covariate_table_varnames <- strsplit(covariate_table_varnames, ",")[[1]]

  # Simplify the data frame
  varnames_to_keep <- c("id", period_varname, "treatment_cohort",
                        outcome_varname, cluster_varname, cohort_varname,
                        sample_varname, covariate_table_varnames)
  if (is_valid_varname(empirical_spec$time_index_varname)) {
    varnames_to_keep <- c(varnames_to_keep, empirical_spec$time_index_varname)
  }
  varnames_to_keep <- Filter(nzchar, varnames_to_keep)
  varnames_to_keep <- Filter(function(x) !is.na(x), varnames_to_keep)

  event_study_df <- event_study_df %>%
    dplyr::select(dplyr::all_of(varnames_to_keep)) %>%
    dplyr::mutate(
      !!outcome_varname := ifelse(
        is.na(!!rlang::sym(outcome_varname)), 2, !!rlang::sym(outcome_varname)
      )
    ) %>%
    dplyr::mutate(time_index = !!rlang::sym(period_varname))

  # If there is a sample variable, filter the data
  if (is_valid_varname(sample_varname)) {
    event_study_df <- event_study_df %>%
      dplyr::filter(!!rlang::sym(sample_varname) == 1)
  }
  # If there is a cohort variable, rename it
  if (is_valid_varname(cohort_varname)) {
    event_study_df <- event_study_df %>%
      dplyr::mutate(treatment_cohort = as.numeric(!!rlang::sym(cohort_varname)))
  }
  event_study_df
}

#' Apply Options from Empirical Specification
#'
#' Updates a list of options by inheriting relevant settings (e.g., standard
#' error method) from an empirical specification object.
#'
#' @param opts A list of options to update. If \code{NULL}, a new list is
#'   created.
#' @param empirical_spec A list or object containing option fields such as
#'   \code{se_method}.
#'
#' @return A list of updated options, including fields inherited from
#'   \code{empirical_spec}.
#'
#' @examples
#' empirical_spec <- list(opts_se_method = "boot_multiplier")
#' apply_opts_from_empirical_spec(list(), empirical_spec)
#'
#' @export
apply_opts_from_empirical_spec <- function(opts, empirical_spec) {
  # Extract standard error method from empirical specification
  opts_se_method <- empirical_spec$opts_se_method

  # Initialize options list if NULL
  if (is.null(opts)) {
    opts <- list()
  }

  # Inherit standard error method if present
  if (is_valid_varname(opts_se_method)) {
    opts$se_method <- opts_se_method
  }

  # Return updated options list
  opts
}
