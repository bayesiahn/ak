#' Format a numeric value with a fixed number of decimal places
#'
#' @param x Numeric value to format
#' @param digits Integer number of decimal places to display
#'
#' @return Character string with the formatted number.
#' @keywords internal
format_fixed_digits <- function(x, digits) {
  if (is.na(x)) return("--")
  sprintf(paste0("%.", digits, "f"), as.numeric(x))
}

#' Get the directory path for saving LaTeX tables
#'
#' @param model List containing model information
#' @param table_type Character string identifying the table type
#' @param create_dir Logical indicating whether to create the directory (default: TRUE)
#'
#' @return Character string with the file path for the LaTeX table.
#' @keywords internal
get_table_file_path <- function(model, table_type, create_dir = TRUE) {
  empirical_spec <- model$empirical_spec
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name

  # Check if the outcome variable name is valid
  J <- length(model$priors)
  dir_name <- paste0("tables/", outcome_varname, "/",
                     outcome_varname, "-", table_type, "-J", J,
                     empirical_spec_suffix(empirical_spec, include_se_info = FALSE), ".tex")

  if (create_dir && !dir.exists(dirname(dir_name))) {
    dir.create(dirname(dir_name), recursive = TRUE)
  }

  return(dir_name)
}

#' Create a LaTeX group header with percentage
#'
#' @param j Integer group index
#' @param pi_j Numeric proportion for the group
#' @param cols Integer number of columns to span (default: 1)
#'
#' @return Character string with the formatted LaTeX header.
#' @keywords internal
format_group_header <- function(j, pi_j, cols = 1) {
  header <- paste0("$j = ", j, "$")
  percentage <- paste0("(", sprintf("%.1f", pi_j*100), "\\%)")

  if (cols > 1) {
    header <- paste0("\\multicolumn{", cols, "}{c}{", header, "}")
    percentage <- paste0("\\multicolumn{", cols, "}{c}{", percentage, "}")
  }

  list(header = header, percentage = percentage)
}

#' Write a LaTeX table to a file
#'
#' @param table Character string containing the LaTeX table code
#' @param file_path Character string with the file path
#'
#' @return Invisibly returns the table string.
#' @keywords internal
write_latex_table <- function(table, file_path) {
  writeLines(table, file_path)
  invisible(table)
}

#' Compute Weighted Mean for a Latent Group
#'
#' Computes the weighted mean of a variable within a specific latent group
#' using the group-specific weights.
#'
#' @param varname Character string specifying the variable name to average.
#' @param unit_df_j Data frame containing unit-level data with a weight_j column.
#'
#' @return Numeric scalar of the weighted mean.
#' @keywords internal
get_j_weighted_mean <- function(varname, unit_df_j) {
  # Extract values
  values <- unit_df_j %>% pull(varname)

  # Extract weights
  weights <- unit_df_j %>% pull(weight_j)

  # Compute group-specific average
  weighted.mean(values, weights, na.rm = TRUE)
}


#' Format summary statistics with appropriate rounding and decimals
#'
#' @param summary_stats Matrix or data frame containing the summary statistics
#' @param digits Integer number of decimal places to round to (default: 2)
#'
#' @return Data frame with formatted numeric values
#'
#' @examples
#' # Format summary statistics with 2 decimal places
#' formatted_stats <- format_summary_stats(stats_matrix)
#' @keywords internal
format_summary_stats <- function(summary_stats, digits = 2) {
  # Convert to data frame if it's a matrix
  summary_df <- as.data.frame(summary_stats)

  # Format all numeric columns to have specified decimal places
  summary_df %>%
    mutate_if(is.numeric, function(x) format(round(x, digits), nsmall = digits))
}

#' Combine summary tables from multiple groups
#'
#' @param model List containing model results and empirical specifications
#' @param tables List of summary statistic matrices, one per group
#'
#' @return Data frame containing the combined summary statistics
#'
#' @examples
#' # Combine summary tables from all groups
#' combined_table <- combine_summary_tables(model, group_tables)
#' @keywords internal
combine_summary_tables <- function(model, tables) {
  # Combine columns from all tables
  covariate_table_names <- model$empirical_spec$covariate_table_names

  # Convert the list of matrices to a single data frame
  combined <- do.call(cbind, tables) %>% as.data.frame()

  # Set row names to variable names
  rownames(combined) <- strsplit(covariate_table_names, ",")[[1]]

  return(combined)
}

#' Convert summary table to long format for LaTeX table generation
#'
#' This function restructures a wide-format summary table with means and standard errors
#' into a long format suitable for creating LaTeX tables.
#'
#' @param data Data frame containing summary table in wide format
#'
#' @return Data frame in long format with Variable, Treatment, Group, Statistic, and Value columns
#'
#' @examples
#' # Convert summary table to long format
#' long_data <- make_summary_table_contents_long(summary_table)
#' @keywords internal
make_summary_table_contents_long <- function(data) {
  data %>%
    tibble::rownames_to_column("Variable") %>%
    tidyr::pivot_longer(-Variable,
                        names_to = c("Treatment", "Statistic", "Group"),
                        names_pattern = "(.*)_(mean|se)_(j\\d+)",
                        values_to = "Value")
}

#' Generate LaTeX column header for group number
#'
#' @param j Integer specifying the group index
#'
#' @return Character string with formatted LaTeX for group number.
#' @keywords internal
get_table_colhead <- function(j) {
  paste0("$J = ", j, "$")
}

#' Generate LaTeX percentage text for group
#'
#' @param pi_j Numeric proportion for the group (between 0 and 1)
#'
#' @return Character string with formatted LaTeX for group percentage.
#' @keywords internal
get_table_percent <- function(pi_j) {
  pi_j_in_percent <- paste0("(", sprintf("%.1f", pi_j*100), "\\%)")
  pi_j_in_percent
}
