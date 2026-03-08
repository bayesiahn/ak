#' Perform weighted t-test comparing variable values between two groups
#'
#' This function conducts a weighted t-test to compare the values of a specified variable
#' between two different groups from the model.
#'
#' @param varname Character string specifying the variable name
#' @param unit_df Data frame containing unit-level data
#' @param model List containing model results including weights
#' @param j_reference Integer specifying the reference group index
#' @param j_compared Integer specifying the comparison group index
#'
#' @return List containing t-statistic, p-value, standard error, and mean difference
#' @export
#'
#' @examples
#' # Compare age between groups 1 and 2
#' age_test <- get_weighted_t_test("age", unit_data, model, 1, 2)
get_weighted_t_test <- function(varname, unit_df, model, j_reference, j_compared) {
  # Extract values and weights for reference group
  weights_mat <- as.matrix(model$weights) %>% hard_classification()

  # Reference group data
  ref_weights <- weights_mat[, j_reference]
  ref_df <- unit_df %>%
    mutate(weight_j = ref_weights[match(id, unique(unit_df$id))]) %>%
    filter(weight_j > 0)
  ref_values <- ref_df %>% pull(varname)
  ref_weights <- ref_df %>% pull(weight_j)
  ref_mean <- weighted.mean(ref_values, ref_weights, na.rm = TRUE)

  # Comparison group data
  comp_weights <- weights_mat[, j_compared]
  comp_df <- unit_df %>%
    mutate(weight_j = comp_weights[match(id, unique(unit_df$id))]) %>%
    filter(weight_j > 0)
  comp_values <- comp_df %>% pull(varname)
  comp_weights <- comp_df %>% pull(weight_j)
  comp_mean <- weighted.mean(comp_values, comp_weights, na.rm = TRUE)

  # Calculate weighted variances
  ref_var <- sum(ref_weights * (ref_values - ref_mean)^2, na.rm = TRUE) /
    (sum(ref_weights, na.rm = TRUE) - sum(ref_weights^2, na.rm = TRUE)/sum(ref_weights, na.rm = TRUE))

  comp_var <- sum(comp_weights * (comp_values - comp_mean)^2, na.rm = TRUE) /
    (sum(comp_weights, na.rm = TRUE) - sum(comp_weights^2, na.rm = TRUE)/sum(comp_weights, na.rm = TRUE))

  # Calculate effective sample sizes
  n_eff_ref <- sum(ref_weights, na.rm = TRUE)^2 / sum(ref_weights^2, na.rm = TRUE)
  n_eff_comp <- sum(comp_weights, na.rm = TRUE)^2 / sum(comp_weights^2, na.rm = TRUE)

  # Standard error of the difference
  se_diff <- sqrt(ref_var/n_eff_ref + comp_var/n_eff_comp)

  # Mean difference
  mean_diff <- ref_mean - comp_mean

  # t-statistic
  t_stat <- mean_diff / se_diff

  # Calculate degrees of freedom (Welch-Satterthwaite approximation)
  df_num <- (ref_var/n_eff_ref + comp_var/n_eff_comp)^2
  df_denom <- (ref_var/n_eff_ref)^2/(n_eff_ref-1) + (comp_var/n_eff_comp)^2/(n_eff_comp-1)
  df <- df_num / df_denom

  # p-value
  p_value <- 2 * pt(-abs(t_stat), df)

  list(
    t_stat = t_stat,
    p_value = p_value,
    se_diff = se_diff,
    mean_diff = mean_diff,
    ref_mean = ref_mean,
    comp_mean = comp_mean,
    ref_se = sqrt(ref_var/n_eff_ref),
    comp_se = sqrt(comp_var/n_eff_comp)
  )
}

#' Generate a table comparing statistics across groups
#'
#' This function creates a data frame containing means, standard errors, and
#' t-test statistics comparing variables across different groups.
#'
#' @param model List containing model results and empirical specifications
#' @param event_study_df Data frame containing event study data
#' @param digits Integer number of decimal places to round to (default: 2)
#' @param perform_hard_classification Logical indicating whether to perform hard classification on weights (default: FALSE)
#'
#' @return Data frame containing cross-group comparison statistics
#' @export
#'
#' @examples
#' # Generate table comparing groups
#' comparison_table <- generate_summary_table_across_types(model, event_study_data)
generate_summary_table_across_types <- function(model, event_study_df, digits = 2,
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
  if (J == 1) {
    # If only one group, just report means and SEs
    results <- data.frame(
      variable = var_names,
      stringsAsFactors = FALSE
    )

    # Compute weights and perform hard classification if needed
    weights <- as.matrix(model$weights)
    if (perform_hard_classification) {
      weights <- weights %>% hard_classification()
    }

    weight_df <- data.frame(id = unique(unit_df$id), weight_j = weights[,1])
    unit_df_weighted <- unit_df %>% left_join(weight_df, by = "id")

    # Calculate means and SEs for each variable
    for (i in seq_along(covariate_table_varnames)) {
      var <- covariate_table_varnames[i]
      results$mean_j1[i] <- get_j_weighted_mean(var, unit_df_weighted)
      results$se_j1[i] <- get_j_weighted_se(var, unit_df_weighted)
    }
  } else {
    # For multiple groups, include pairwise comparisons with group 1 as reference
    results <- data.frame(
      variable = var_names,
      stringsAsFactors = FALSE
    )

    # Calculate statistics for each group
    for (j in 1:J) {
      # Compute weights for this group
      weights <- as.matrix(model$weights) %>% hard_classification()
      weight_df <- data.frame(id = unique(unit_df$id), weight_j = weights[,j])
      unit_df_j <- unit_df %>% left_join(weight_df, by = "id")

      # Calculate means and SEs for each variable
      for (i in seq_along(covariate_table_varnames)) {
        var <- covariate_table_varnames[i]
        col_mean <- paste0("mean_j", j)
        col_se <- paste0("se_j", j)
        results[[col_mean]][i] <- get_j_weighted_mean(var, unit_df_j)
        results[[col_se]][i] <- get_j_weighted_se(var, unit_df_j)
      }
    }

    # Calculate difference statistics for each variable (compare all groups to group 1)
    for (j in 2:J) {
      for (i in seq_along(covariate_table_varnames)) {
        var <- covariate_table_varnames[i]
        test_result <- get_weighted_t_test(var, unit_df, model, 1, j)

        col_diff <- paste0("diff_j1_j", j)
        col_diff_se <- paste0("diff_se_j1_j", j)
        col_t_stat <- paste0("t_stat_j1_j", j)

        results[[col_diff]][i] <- test_result$mean_diff
        results[[col_diff_se]][i] <- test_result$se_diff
        results[[col_t_stat]][i] <- test_result$t_stat
      }
    }
  }

  # Format numeric columns
  results_formatted <- results %>%
    mutate_if(is.numeric, function(x) round(x, digits = digits))

  return(results_formatted)
}

#' Generate a LaTeX table comparing statistics across groups
#'
#' This function creates a formatted LaTeX table showing means, standard errors,
#' and difference statistics comparing variables across different groups,
#' with all numeric values displayed with a fixed number of decimal places.
#'
#' @param model List containing model results and empirical specifications
#' @param event_study_df Data frame containing event study data
#' @param digits Integer number of decimal places to display (default: 2)
#' @param diff_col_name Character string for the difference column header (default: "Difference")
#'
#' @return Character string containing the LaTeX table code
#' @export
generate_summary_table_across_types_tex <- function(model, event_study_df,
                                                    digits = 2,
                                                    diff_col_name = "Difference") {
  # Get cross-group comparison statistics
  results <- generate_summary_table_across_types(model, event_study_df, digits)

  # Number of groups
  J <- length(model$priors)

  # Start building LaTeX table
  header <- "\\centering\n"

  if (J == 1) {
    # For single group, just show means and SEs
    header <- paste0(header, "\\begin{tabular}{lc}\n\\hline\\hline\n")

    group_format <- format_group_header(1, model$priors[1])
    header <- paste0(header, " Transition Type (Share) & ", group_format$header, " \\\\\n")
    header <- paste0(header, " & ", group_format$percentage, " \\\\\n\\hline\n")

    # Body
    body <- ""
    for (i in 1:nrow(results)) {
      var_name <- results$variable[i]

      # Format with fixed decimal places
      mean_j1 <- format_fixed_digits(results$mean_j1[i], digits)
      se_j1 <- format_fixed_digits(results$se_j1[i], digits)

      # Add variable row
      body <- paste0(body, var_name, " & ", mean_j1, " \\\\\n")

      # Add SE row
      body <- paste0(body, " & (", se_j1, ") \\\\\n")
    }
  } else {
    # For multiple groups with comparisons
    # Determine number of columns based on groups
    num_cols <- J + 1  # J group columns + 1 difference column

    # Create tabular environment
    header <- paste0(header, "\\begin{tabular}{l", paste(rep("c", num_cols), collapse = ""), "}\n")
    header <- paste0(header, "\\hline\\hline\n")

    # Column headers for groups and difference
    group_headers <- sapply(1:J, function(j) {
      format_group_header(j, model$priors[j])$header
    })

    header <- paste0(header, " Transition Type (Share) & ",
                     paste(c(group_headers, diff_col_name), collapse = " & "), " \\\\\n")

    # Group percentages row
    percentages <- sapply(1:J, function(j) {
      format_group_header(j, model$priors[j])$percentage
    })

    header <- paste0(header, " & ", paste(c(percentages, ""), collapse = " & "), " \\\\\n\\hline\n")

    # Body
    body <- ""
    for (i in 1:nrow(results)) {
      var_name <- results$variable[i]

      # Create row with means and difference
      means <- sapply(1:J, function(j) {
        format_fixed_digits(results[[paste0("mean_j", j)]][i], digits)
      })

      # If more than one group, add difference column
      if (J > 1) {
        t_stat <- format_fixed_digits(results[[paste0("diff_j1_j", 2)]][i], digits)
        row_values <- c(means, t_stat)
      } else {
        row_values <- means
      }

      # Add variable row with means and t-statistic
      body <- paste0(body, var_name, " & ", paste(row_values, collapse = " & "), " \\\\\n")

      # Create row with standard errors
      ses <- sapply(1:J, function(j) {
        paste0("(", format_fixed_digits(results[[paste0("se_j", j)]][i], digits), ")")
      })

      # If more than one group, add difference SE
      if (J > 1) {
        diff_se <- format_fixed_digits(results[[paste0("diff_se_j1_j", 2)]][i], digits)
        se_values <- c(ses, paste0("(", diff_se, ")"))
      } else {
        se_values <- ses
      }

      # Add SE row
      body <- paste0(body, " & ", paste(se_values, collapse = " & "), " \\\\\n")
    }
  }

  # Complete the table
  table <- paste0(header, body, "\\hline\\hline\n\\end{tabular}\n")

  # Save table to file
  file_path <- get_table_file_path(model, "across-types")
  write_latex_table(table, file_path)

  return(table)
}
