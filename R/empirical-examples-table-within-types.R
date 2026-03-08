#' Perform weighted t-test comparing treated and control units within a group
#'
#' This function conducts a weighted t-test to compare treated and control units
#' for a specified variable within a single group.
#'
#' @param varname Character string specifying the variable name
#' @param unit_df Data frame containing unit-level data with weight_j column
#' @param treatment_var Character string specifying the treatment indicator variable (default: "treatment_cohort")
#'
#' @return List containing t-statistic, p-value, standard error, and mean difference
#' @export
#'
#' @examples
#' # Compare income between treated and control units in group 1
#' income_test <- get_within_type_t_test("income", group1_df)
get_within_type_t_test <- function(varname, unit_df, treatment_var = "treatment_cohort") {
  # Split by treatment status
  treated_df <- unit_df %>% filter(!!sym(treatment_var) != 0)
  control_df <- unit_df %>% filter(!!sym(treatment_var) == 0)

  # Get values and weights
  treated_values <- treated_df %>% pull(varname)
  treated_weights <- treated_df %>% pull(weight_j)
  control_values <- control_df %>% pull(varname)
  control_weights <- control_df %>% pull(weight_j)

  # Calculate weighted means
  treated_mean <- weighted.mean(treated_values, treated_weights, na.rm = TRUE)
  control_mean <- weighted.mean(control_values, control_weights, na.rm = TRUE)

  # Calculate weighted variances
  treated_var <- sum(treated_weights * (treated_values - treated_mean)^2, na.rm = TRUE) /
    (sum(treated_weights, na.rm = TRUE) - sum(treated_weights^2, na.rm = TRUE)/sum(treated_weights, na.rm = TRUE))

  control_var <- sum(control_weights * (control_values - control_mean)^2, na.rm = TRUE) /
    (sum(control_weights, na.rm = TRUE) - sum(control_weights^2, na.rm = TRUE)/sum(control_weights, na.rm = TRUE))

  # Calculate effective sample sizes
  n_eff_treated <- sum(treated_weights, na.rm = TRUE)^2 / sum(treated_weights^2, na.rm = TRUE)
  n_eff_control <- sum(control_weights, na.rm = TRUE)^2 / sum(control_weights^2, na.rm = TRUE)

  # Standard error of the difference
  se_diff <- sqrt(treated_var/n_eff_treated + control_var/n_eff_control)

  # Mean difference (treated - control)
  mean_diff <- treated_mean - control_mean

  # t-statistic
  t_stat <- mean_diff / se_diff

  # Calculate degrees of freedom (Welch-Satterthwaite approximation)
  df_num <- (treated_var/n_eff_treated + control_var/n_eff_control)^2
  df_denom <- (treated_var/n_eff_treated)^2/(n_eff_treated-1) + (control_var/n_eff_control)^2/(n_eff_control-1)
  df <- df_num / df_denom

  # p-value
  p_value <- 2 * pt(-abs(t_stat), df)

  list(
    t_stat = t_stat,
    p_value = p_value,
    se_diff = se_diff,
    mean_diff = mean_diff,
    treated_mean = treated_mean,
    control_mean = control_mean,
    treated_se = sqrt(treated_var/n_eff_treated),
    control_se = sqrt(control_var/n_eff_control)
  )
}

#' Generate a table of treatment effects within each group
#'
#' This function creates a data frame showing treatment effects and statistics
#' for treated versus control units within each group.
#'
#' @param model List containing model results and empirical specifications
#' @param event_study_df Data frame containing event study data
#' @param treatment_var Character string specifying the treatment indicator variable (default: "treatment_cohort")
#' @param digits Integer number of decimal places to round to (default: 2)
#' @param perform_hard_classification Logical indicating whether to perform hard classification on weights (default: FALSE)
#'
#' @return Data frame containing within-group treatment effect statistics
#' @export
#'
#' @examples
#' # Generate table of within-group treatment effects
#' treatment_effects <- generate_summary_table_within_type(model, event_study_data)
generate_summary_table_within_type <- function(model, event_study_df,
                                                treatment_var = "treatment_cohort",
                                                digits = 2,
                                               perform_hard_classification = FALSE) {
  # Get empirical specifications
  empirical_spec <- model$empirical_spec
  covariate_table_varnames <- empirical_spec$covariate_table_varnames
  covariate_table_varnames <- strsplit(covariate_table_varnames, ",")[[1]]
  covariate_table_names <- empirical_spec$covariate_table_names
  var_names <- strsplit(covariate_table_names, ",")[[1]]
  reference_time_index <- empirical_spec$covariate_table_reference_time_index

  # Prepare data
  unit_df <- event_study_df %>%
    apply_sample_restriction(empirical_spec) %>%
    filter(time_index == reference_time_index)

  # Number of groups
  J <- length(model$priors)

  # Initialize results data frame
  results <- data.frame(
    variable = var_names,
    stringsAsFactors = FALSE
  )

  # For each group, calculate treatment effects
  for (j in 1:J) {
    # Compute weights for this group
    weights <- as.matrix(model$weights)
    if (perform_hard_classification) {
      weights <- weights %>% hard_classification()
    }
    weight_df <- data.frame(id = unique(unit_df$id), weight_j = weights[,j])
    unit_df_j <- unit_df %>% left_join(weight_df, by = "id") %>% filter(weight_j > 0)

    # For each variable, calculate treatment effect statistics
    for (i in seq_along(covariate_table_varnames)) {
      var <- covariate_table_varnames[i]
      test_result <- get_within_type_t_test(var, unit_df_j, treatment_var)

      # Store results
      col_prefix <- paste0("j", j, "_")
      results[[paste0(col_prefix, "treated_mean")]][i] <- test_result$treated_mean
      results[[paste0(col_prefix, "treated_se")]][i] <- test_result$treated_se
      results[[paste0(col_prefix, "control_mean")]][i] <- test_result$control_mean
      results[[paste0(col_prefix, "control_se")]][i] <- test_result$control_se
      results[[paste0(col_prefix, "diff")]][i] <- test_result$mean_diff
      results[[paste0(col_prefix, "diff_se")]][i] <- test_result$se_diff
      results[[paste0(col_prefix, "t_stat")]][i] <- test_result$t_stat
      results[[paste0(col_prefix, "p_value")]][i] <- test_result$p_value
    }
  }

  # Format numeric columns
  results_formatted <- results %>%
    mutate_if(is.numeric, function(x) round(x, digits = digits))

  return(results_formatted)
}

#' Generate a LaTeX table of treatment effects within each group
#'
#' This function creates a formatted LaTeX table showing treatment effects
#' for treated versus control units within each group, with all numeric values
#' displayed with a fixed number of decimal places.
#'
#' @param model List containing model results and empirical specifications
#' @param event_study_df Data frame containing event study data
#' @param treatment_var Character string specifying the treatment indicator variable (default: "treatment_cohort")
#' @param digits Integer number of decimal places to display (default: 2)
#' @param diff_col_name Character string for the difference column header (default: "Diff.")
#'
#' @return Character string containing the LaTeX table code
#' @export
generate_summary_table_within_type_tex <- function(model, event_study_df,
                                                    treatment_var = "treatment_cohort",
                                                    digits = 2,
                                                    diff_col_name = "Difference") {
  # Get within-group treatment effect statistics
  results <- generate_summary_table_within_type(model, event_study_df, treatment_var, digits)

  # Number of groups
  J <- length(model$priors)

  # Start building LaTeX table
  header <- "\\centering\n"

  # Determine number of columns for tabular environment (3 per group: Treated, Control, Diff.)
  col_spec <- paste0("l", paste(rep("ccc", J), collapse = ""))
  header <- paste0(header, "\\begin{tabular}{", col_spec, "}\n")
  header <- paste0(header, "\\hline\\hline\n")

  # Column headers for transition type and groups
  group_col_headers <- c()
  group_percentages <- c()

  for (j in 1:J) {
    group_format <- format_group_header(j, model$priors[j], cols = 3)
    group_col_headers <- c(group_col_headers, group_format$header)
    group_percentages <- c(group_percentages, group_format$percentage)
  }

  header <- paste0(header, " Transition Type (Share) & ",
    paste(group_col_headers, collapse = " & "), " \\\\\n")
  header <- paste0(header, " & ",
    paste(group_percentages, collapse = " & "), " \\\\\n")

  # Treatment status headers
  treatment_headers <- c()
  for (j in 1:J) {
    treatment_header <- paste0("Treated & Control & ", diff_col_name)
    treatment_headers <- c(treatment_headers, treatment_header)
  }

  header <- paste0(header, " & ", paste(treatment_headers, collapse = " & "), " \\\\\n\\hline\n")

  # Body
  body <- ""
  for (i in 1:nrow(results)) {
    var_name <- results$variable[i]

    # First row: Variable name and means/differences
    means_row <- c()
    for (j in 1:J) {
      col_prefix <- paste0("j", j, "_")

      # Format values with fixed decimal places
      treated_mean <- format_fixed_digits(results[[paste0(col_prefix, "treated_mean")]][i], digits)
      control_mean <- format_fixed_digits(results[[paste0(col_prefix, "control_mean")]][i], digits)
      diff <- format_fixed_digits(results[[paste0(col_prefix, "diff")]][i], digits)

      group_means <- paste0(treated_mean, " & ", control_mean, " & ", diff)
      means_row <- c(means_row, group_means)
    }

    # Add variable row with means and differences
    body <- paste0(body, var_name, " & ", paste(means_row, collapse = " & "), " \\\\\n")

    # Second row: Blank and standard errors
    se_row <- c()
    for (j in 1:J) {
      col_prefix <- paste0("j", j, "_")

      # Format standard errors with fixed decimal places
      treated_se <- format_fixed_digits(results[[paste0(col_prefix, "treated_se")]][i], digits)
      control_se <- format_fixed_digits(results[[paste0(col_prefix, "control_se")]][i], digits)
      diff_se <- format_fixed_digits(results[[paste0(col_prefix, "diff_se")]][i], digits)

      group_ses <- paste0("(", treated_se, ") & (", control_se, ") & (", diff_se, ")")
      se_row <- c(se_row, group_ses)
    }

    # Add SE row
    body <- paste0(body, " & ", paste(se_row, collapse = " & "), " \\\\\n")
  }

  # Complete the table
  table <- paste0(header, body, "\\hline\\hline\n\\end{tabular}\n")

  # Save table to file
  file_path <- get_table_file_path(model, "within-types")
  write_latex_table(table, file_path)

  return(table)
}
