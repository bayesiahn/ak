#' Calculate weighted standard error of a variable
#'
#' This function computes the weighted standard error for a variable using
#' the appropriate formula for weighted data.
#'
#' @param varname Character string specifying the variable name
#' @param unit_df_j Data frame containing unit-level data with weight_j column
#'
#' @return Numeric weighted standard error of the variable
#' @export
#'
#' @examples
#' # Calculate weighted standard error of income for group 1
#' income_se_group1 <- get_j_weighted_se("income", group1_data)
get_j_weighted_se <- function(varname, unit_df_j) {
  # Extract values
  values <- unit_df_j %>% pull(varname)
  # Extract weights
  weights <- unit_df_j %>% pull(weight_j)

  # Return NA if insufficient data
  if (length(values) <= 1 || all(is.na(values))) {
    return(NA)
  }

  # Calculate weighted mean
  w_mean <- weighted.mean(values, weights, na.rm = TRUE)

  # Calculate weighted variance
  # Formula: sum(w_i * (x_i - w_mean)^2) / (sum(w_i) - sum(w_i^2)/sum(w_i))
  numerator <- sum(weights * (values - w_mean)^2, na.rm = TRUE)
  denominator <- sum(weights, na.rm = TRUE) - sum(weights^2, na.rm = TRUE)/sum(weights, na.rm = TRUE)

  if (denominator <= 0) {
    return(NA)
  }

  w_var <- numerator / denominator

  # Standard error = sqrt(weighted variance / n_eff)
  # n_eff = (sum(w_i))^2 / sum(w_i^2)
  n_eff <- sum(weights, na.rm = TRUE)^2 / sum(weights^2, na.rm = TRUE)

  w_se <- sqrt(w_var / n_eff)
  return(w_se)
}

#' Get summary statistics for a specific group
#'
#' This function calculates summary statistics for covariates in a specific group,
#' separating results for treated and control units.
#'
#' @param j Integer specifying the group index
#' @param model List containing model results and empirical specifications
#' @param event_study_df Data frame containing event study data
#' @param perform_hard_classification Logical indicating whether to perform hard classification on weights (default: FALSE)
#'
#' @return Matrix of summary statistics for treated, control, and all units in group j
#' @export
#'
#' @examples
#' # Get summary statistics for group 1
#' group1_stats <- generate_summary_table_j(1, model, event_study_data)
generate_summary_table_j <- function(j, model, event_study_df,
                                     perform_hard_classification = FALSE) {
  # Retrieve covariates to be taken
  empirical_spec <- model$empirical_spec
  covariate_table_varnames <- empirical_spec$covariate_table_varnames
  covariate_table_varnames <- strsplit(covariate_table_varnames, ",")[[1]]
  reference_time_index <- empirical_spec$covariate_table_reference_time_index

  # Restrict the sample to the time index of interest
  unit_df_j <- event_study_df %>%
    apply_sample_restriction(empirical_spec) %>%
    filter(time_index == reference_time_index)

  # Compute weights and perform hard classification if needed
  weights <- as.matrix(model$weights)

  if (perform_hard_classification) {
    weights <- weights %>% hard_classification()
  }
  weight_df_j <- data.frame(id = unique(unit_df_j$id),
                            weight_j = weights[,j])
  unit_df_j <- unit_df_j %>%
    left_join(weight_df_j, by = "id")

  # Split by treatment status
  unit_df_j_treated <- unit_df_j %>% filter(treatment_cohort != 0)
  unit_df_j_control <- unit_df_j %>% filter(treatment_cohort == 0)

  # Calculate means
  treated_means <- sapply(covariate_table_varnames, get_j_weighted_mean, unit_df_j_treated)
  control_means <- sapply(covariate_table_varnames, get_j_weighted_mean, unit_df_j_control)
  all_means <- sapply(covariate_table_varnames, get_j_weighted_mean, unit_df_j)

  # Calculate standard errors
  treated_se <- sapply(covariate_table_varnames, get_j_weighted_se, unit_df_j_treated)
  control_se <- sapply(covariate_table_varnames, get_j_weighted_se, unit_df_j_control)
  all_se <- sapply(covariate_table_varnames, get_j_weighted_se, unit_df_j)

  # Combine means and SEs in a matrix
  # Format: each variable has 6 columns - 3 means and 3 SEs
  result_matrix <- matrix(NA, nrow = length(covariate_table_varnames), ncol = 6)
  result_matrix[, 1] <- treated_means
  result_matrix[, 2] <- treated_se
  result_matrix[, 3] <- control_means
  result_matrix[, 4] <- control_se
  result_matrix[, 5] <- all_means
  result_matrix[, 6] <- all_se

  # Set column names
  colnames(result_matrix) <- c(
    paste0("treated_mean_j", j),
    paste0("treated_se_j", j),
    paste0("control_mean_j", j),
    paste0("control_se_j", j),
    paste0("all_mean_j", j),
    paste0("all_se_j", j)
  )

  # Set row names
  rownames(result_matrix) <- covariate_table_varnames

  return(result_matrix)
}

#' Generate a complete summary table for all groups
#'
#' This function combines summary statistics for all groups into a single data frame,
#' formatting numeric values to show two decimal places.
#'
#' @param model List containing model results and empirical specifications
#' @param event_study_df Data frame containing event study data
#' @param digits Integer number of decimal places to round to (default: 2)
#'
#' @return Data frame containing summary statistics for all groups, formatted with specified decimal places
#' @export
#'
#' @examples
#' # Generate complete summary table
#' full_summary <- generate_summary_table(model, event_study_data)
generate_summary_table <- function(model, event_study_df, digits = 2,
                                   perform_hard_classification = FALSE) {
  # Calculate summary statistics for each group
  tables <- lapply(1:length(model$priors), generate_summary_table_j, model, event_study_df,
                   perform_hard_classification = perform_hard_classification)

  # Combine tables from all groups
  combined <- combine_summary_tables(model, tables)

  # Format summary statistics
  formatted <- format_summary_stats(combined, digits)

  return(formatted)
}

#' Generate a LaTeX table of summary statistics by group and treatment status
#'
#' This function creates a formatted LaTeX table showing summary statistics for
#' each covariate by group and treatment status. The table includes group proportions
#' and separates treated and control units, with standard errors shown in a separate row
#' below each mean value.
#'
#' @param model List containing model results and empirical specifications
#' @param event_study_df Data frame containing event study data
#' @param var_label Character string specifying the variable column name (default: "Variable")
#' @param group_label Character string specifying the group column name (default: "Group")
#' @param treatment_label Character string specifying the treatment column name (default: "Treatment")
#' @param value_label Character string specifying the value column name (default: "Value")
#' @param statistic_label Character string specifying the statistic column name (default: "Statistic")
#'
#' @return Character string containing the LaTeX table code
#' @export
generate_summary_table_tex <- function(model, event_study_df,
                                       var_label = "Variable",
                                       group_label = "Group",
                                       treatment_label = "Treatment",
                                       value_label = "Value",
                                       statistic_label = "Statistic") {

  # Generate summary table and convert to long format
  data <- generate_summary_table(model, event_study_df) %>%
    make_summary_table_contents_long()

  groups <- unique(data[[group_label]])
  treatments <- unique(data[[treatment_label]])

  # Header
  J <- length(model$priors)
  j_in_models <- 1:J

  header <- paste0("\n\\centering\n\\begin{tabular}{l",
                   paste(rep("cc", length(groups)), collapse = ""),
                   "}\n\\hline\\hline\n")

  # First row with group numbers
  header <- paste0(header, " Transition Type (Share) & ",
                   paste(sapply(1:J, function(j) {
                     paste0("\\multicolumn{2}{c}{", get_table_colhead(j), "}")
                   }), collapse = " & "), " \\\\\n")

  # Second row with percentages
  header <- paste0(header, " & ",
                   paste(sapply(1:J, function(j) {
                     paste0("\\multicolumn{2}{c}{", get_table_percent(model$priors[j_in_models[j]]), "}")
                   }), collapse = " & "), " \\\\\n")

  # Treatment status row
  sub_header <- " Treatment Status & "
  sub_header <- paste0(sub_header,
                       paste(rep("\\multicolumn{1}{c}{Treated} & \\multicolumn{1}{c}{Control}",
                                 length(groups)),
                             collapse = " & "),
                       " \\\\\n\\hline\n")

  # Body
  body <- ""
  for (var in unique(data[[var_label]])) {
    var_data <- data[data[[var_label]] == var, ]

    # First row: Variable name and mean values
    means_row <- character()
    se_row <- character()

    for (g in groups) {
      group_data <- var_data[var_data[[group_label]] == g, ]

      # Get treated mean and SE
      treated_mean <- group_data[group_data[[treatment_label]] == "treated" &
                                   group_data[[statistic_label]] == "mean", value_label]
      treated_se <- group_data[group_data[[treatment_label]] == "treated" &
                                 group_data[[statistic_label]] == "se", value_label]

      # Get control mean and SE
      control_mean <- group_data[group_data[[treatment_label]] == "control" &
                                   group_data[[statistic_label]] == "mean", value_label]
      control_se <- group_data[group_data[[treatment_label]] == "control" &
                                 group_data[[statistic_label]] == "se", value_label]

      # Add to means row and SE row
      means_row <- c(means_row, treated_mean, control_mean)
      se_row <- c(se_row, paste0("(", treated_se, ")"), paste0("(", control_se, ")"))
    }

    # Add mean row to body
    means_text <- paste(c(var, means_row), collapse = " & ")
    body <- paste0(body, means_text, " \\\\\n")

    # Add SE row to body (blank in first column)
    se_text <- paste(c("", se_row), collapse = " & ")
    body <- paste0(body, se_text, " \\\\\n")
  }

  # Combine
  table <- paste0(header, sub_header, body, "\\hline\\hline\n\\end{tabular}\n")

  # Save in get_summary_table_directory
  dir_name <- get_table_file_path(model, "within-types-without-diffs")
  if (!dir.exists(dirname(dir_name))) {
    dir.create(dirname(dir_name), recursive = TRUE)
  }
  writeLines(table, dir_name)

  table
}
