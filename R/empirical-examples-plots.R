PLOT_YLIM_AVERAGES_DIFF <- c(0, 1)
PLOT_YLIM_AVERAGES <- c(0, 1)

#' Get Y values in numeric and G from Model
#'
#' Extracts outcome (`y`) as numeric values and treatment cohort (`g`) data from the model's event study data frame.
#'
#' @param model A model object containing the event study data.
#' @return A list containing `y_numeric` and `g` variables extracted from the data.
#' @export
get_y_numeric_and_g_from_model <- function(model) {
  y_g <- long_df_to_y_and_g(
    model$event_study_df,
    model$empirical_spec$outcome_varname,
    tname = "time_index",
    idname = "id",
    gname = "treatment_cohort"
  )

  # Convert treatment cohort to 0-indexed
  y_g$y_numeric <- matrix(model$Y$values[factor(y_g$y)], ncol = ncol(y_g$y))

  return(y_g)
}

#' Plot Averages
#'
#' Generates plots for average outcomes, optionally by treatment cohort difference.
#'
#' @param model A model object containing the event study data.
#' @param difference_by_treatment A logical value. If TRUE, shows differences
#'   by treatment cohort.
#' @param save A logical value. If TRUE (default), saves the plot to file.
#' @return A ggplot object containing the average plots.
#' @export
plot_averages <- function(model,
                          difference_by_treatment = FALSE,
                          save = TRUE) {
  # Retrieve outcome name and variable name
  empirical_spec <- model$empirical_spec
  outcome_name <- empirical_spec$outcome_name
  outcome_varname <- empirical_spec$outcome_varname
  y_lab_average_1 <- empirical_spec$y_lab_average_1
  x_lab <- empirical_spec$x_lab

  # Retrieve y and g from the model
  y_g <- get_y_numeric_and_g_from_model(model)
  y <- y_g$y_numeric
  g <- y_g$g

  p <- plot_transition_averages(
    y, g, model$weights,
    difference_by_treatment = difference_by_treatment, display_ci = FALSE) +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab_average_1)

  p_with_title <- p +
    ggplot2::ggtitle(empirical_spec_to_plot_title(
      empirical_spec, "Type-Specific Trends"
    ))

  # Save plots if requested
  if (save) {
    diff_suffix <- ifelse(difference_by_treatment, "-diff-by-cohort", "")
    save_plot(get_plot_averages_directory(
      empirical_spec, J = length(model$priors), suffix = diff_suffix), p)
    save_plot(get_plot_averages_directory(
      empirical_spec, J = length(model$priors),
      suffix = paste0(diff_suffix, "-with-title")), p_with_title)
  }

  p_with_title
}

#' Plot LTATT
#'
#' Generates a plot for Group Average Treatment Effects (LTATT).
#'
#' @param model A model object containing the event study data.
#' @param display_ci Determines if confidence intervals should be displayed.
#'   Default is FALSE.
#' @param uniform_ci Logical; if TRUE (default), use uniform (simultaneous)
#'   confidence intervals; if FALSE, use pointwise intervals.
#' @param save A logical value. If TRUE (default), saves the plot to file.
#' @return A ggplot object containing the LTATT plot.
#' @export
plot_ltatt_and_save <- function(model, display_ci = FALSE, uniform_ci = TRUE,
                                 save = TRUE) {
  # Retrieve outcome name and variable name
  empirical_spec <- model$empirical_spec
  outcome_name <- empirical_spec$outcome_name
  outcome_varname <- empirical_spec$outcome_varname
  x_lab <- empirical_spec$x_lab

  p_gatt <- plot_ltatt(model, display_ci = display_ci, uniform_ci = uniform_ci) +
    ggplot2::xlab(x_lab)

  if (save) {
    save_plot(get_plot_ltatt_directory(empirical_spec,
                                       J = length(model$priors)), p_gatt)
  }

  p_gatt +
    ggplot2::ggtitle(empirical_spec_to_plot_title(
      empirical_spec, "Type-Specific ATT"
    ))
}

#' Plot LTATT by pretreatment outcomes
#'
#' Generates a plot for Group Average Treatment Effects (LTATT) by
#' pretreatment outcomes (outcomes just before intervention).
#'
#' @param model A model object containing the event study data.
#' @param save A logical value. If TRUE (default), saves the plot to file.
#' @return A ggplot object containing the LTATT plot.
#' @export
plot_ltatt_by_pretreatment_outcomes_and_save <- function(model, save = TRUE) {
  # Retrieve outcome name and variable name
  empirical_spec <- model$empirical_spec
  outcome_name <- empirical_spec$outcome_name
  outcome_varname <- empirical_spec$outcome_varname
  x_lab <- empirical_spec$x_lab

  p_gatt <- plot_ltatt_by_pretreatment_outcomes(model) +
    ggplot2::xlab(x_lab)

  if (save) {
    save_plot(get_plot_ltatt_by_pretreatment_outcomes_directory(
      empirical_spec, J = length(model$priors)), p_gatt)
  }

  p_gatt +
    ggplot2::ggtitle(empirical_spec_to_plot_title(
      empirical_spec, "Type-Specific ATT by Pretreatment Outcomes"
    ))
}

#' Plot Transitions
#'
#' Generates plots for transition probabilities between outcome states.
#'
#' @param model A model object containing the event study data.
#' @param display_ci A logical value. If TRUE, confidence intervals are displayed.
#' @param uniform_ci A logical value. If TRUE, uses uniform confidence intervals.
#' @param difference_by_treatment A logical value. If TRUE, shows differences
#'   by treatment cohort.
#' @param pretreatment_only A logical value. If TRUE, only pretreatment periods
#'   are plotted.
#' @param y_lab_function Optional custom label function for the y-axis.
#'   Takes y_name_from and y_name_to to return ylab. If NULL, defaults are used.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param save A logical value. If TRUE (default), saves the plots to file.
#' @return A list of ggplot objects containing the transition plots.
#' @export
plot_transitions <- function(model,
                             display_ci = TRUE,
                             uniform_ci = TRUE,
                             difference_by_treatment = FALSE,
                             pretreatment_only = FALSE,
                             y_lab_function = NULL,
                             width = 9, height = 9,
                             save = TRUE) {
  # Retrieve outcome name and variable name
  empirical_spec <- model$empirical_spec
  outcome_name <- empirical_spec$outcome_name
  outcome_varname <- empirical_spec$outcome_varname

  # Retrieve x_lab and y_lab from the model
  x_lab <- model$empirical_spec$x_lab

  # Extract variable names
  y_names <- model$Y$names
  plot_title <- empirical_spec_to_plot_title(
    empirical_spec, "Type-Specific Transition Probabilities")
  suffix <- ifelse(difference_by_treatment, "-diff", "")
  suffix <- paste0(suffix, ifelse(pretreatment_only, "-pretreatment", ""))
  suffix <- paste0(suffix, ifelse(uniform_ci, "-uniform_ci", ""))

  ps <- list()
  for (y_name_from in y_names) {
    for (y_name_to in y_names) {
      plot_name <- get_transition_plot_suffix(y_name_from, y_name_to)
      y_lab <- sprintf("%s-to-%s Rates", y_name_from, y_name_to)
      if (!is.null(y_lab_function)) {
        y_lab <- y_lab_function(y_name_from, y_name_to)
      }
      p <- plot_transition_probabilities(
        model, y_name_to = y_name_to, y_name_from = y_name_from,
        display_ci = display_ci, pretreatment_only = pretreatment_only,
        uniform_ci = uniform_ci
      ) +
        ggplot2::xlab(x_lab) +
        ggplot2::ylab(y_lab)
      p_with_title <- p + ggplot2::ggtitle(plot_title)

      # Save plot if requested
      if (save) {
        save_plot(get_plot_transitions_directory(empirical_spec,
          J = length(model$priors), suffix = paste0(suffix, plot_name)), p)
        save_plot(get_plot_transitions_directory(empirical_spec,
          J = length(model$priors), suffix = paste0(suffix, plot_name, "-with-title")),
          p_with_title)
      }

      ps[[plot_name]] <- p
    }
  }

  # Generate representative plots
  plotlist <- .make_transition_plotlist(ps, y_names)
  p <- ggpubr::ggarrange(
    plotlist = plotlist,
    ncol = length(y_names),
    nrow = length(y_names),
    common.legend = TRUE,
    legend = "bottom"
  ) +
    guides(
      nrow = 1,
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2),
      shape = guide_legend(order = 2)
    )

  p_with_title <- p + ggplot2::ggtitle(plot_title)

  # Save plots if requested
  if (save) {
    save_plot(get_plot_transitions_directory(empirical_spec,
      J = length(model$priors), suffix = suffix), p, width = width, height = height)
    save_plot(get_plot_transitions_directory(empirical_spec,
      J = length(model$priors), suffix = paste0(suffix, "-with-title")),
      p_with_title, width = width, height = height)
  }

  ps
}

# helper: row-major plotlist (rows = y_name_to, cols = y_name_from)
.make_transition_plotlist <- function(ps, y_names) {
  unlist(
    lapply(y_names, function(y_to)
      lapply(y_names, function(y_from)
        ps[[get_transition_plot_suffix(y_name_from = y_from, y_name_to = y_to)]]
      )
    ),
    recursive = FALSE
  )
}

get_transition_plot_suffix <- function(y_name_from, y_name_to) {
  paste0("-from-", y_name_from, "-to-", y_name_to)
}

#' Plot Data Trends
#'
#' Generates a plot showing data trends for the outcome variable over time,
#' relative to a policy change.
#'
#' @param model A model object containing the event study data.
#' @param display_ci Logical. Whether to display confidence intervals. Default
#'   is TRUE.
#' @param ylim_min Minimum value for the y-axis. Default is -0.1.
#' @param ylim_max Maximum value for the y-axis. Default is 1.1.
#' @param save A logical value. If TRUE (default), saves the plot to file.
#' @return A ggplot object displaying the data trends.
#' @export
plot_data_trends <- function(model,
                             display_ci = TRUE,
                             ylim_min = -0.1, ylim_max = 1.1,
                             save = TRUE) {
  # Extract outcome name and variable name
  empirical_spec <- model$empirical_spec
  outcome_name <- model$empirical_spec$outcome_name
  outcome_varname <- model$empirical_spec$outcome_varname
  x_lab <- model$empirical_spec$x_lab
  J <- length(model$priors)

  # Extract data
  y_g <- get_y_numeric_and_g_from_model(model)
  y_numeric <- y_g$y_numeric
  g <- y_g$g
  p <- plot_transition_averages(y_numeric, g, difference_by_treatment = FALSE,
                          display_ci = display_ci) +
    ggplot2::scale_color_manual(values = "Black", name = "", guide = "none") +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(outcome_name) +
    ggplot2::ylim(ylim_min, ylim_max)

  p_with_title <- p +
    ggplot2::ggtitle(empirical_spec_to_plot_title(
      empirical_spec, "Trends", display_se_info = FALSE
    ))

  # Save plots if requested
  if (save) {
    save_plot(get_plot_averages_directory(empirical_spec), p)
    save_plot(get_plot_averages_directory(empirical_spec,
      suffix = "-with-title"), p_with_title)
  }

  p_with_title
}


#' Plot Data Trends (All Outcomes)
#'
#' Generates a plot showing data trends for all outcome states over time,
#' relative to a policy change.
#'
#' @param model A model object containing the event study data.
#' @param display_ci Logical. Whether to display confidence intervals. Default
#'   is TRUE.
#' @param recenter_t Logical. Whether to recenter time variable t. Default is
#'   TRUE.
#' @param display_vertical_event_lines Logical. Whether to display vertical
#'   event lines. Default is TRUE.
#' @param ylim_min Minimum value for the y-axis. Default is -0.1.
#' @param ylim_max Maximum value for the y-axis. Default is 1.1.
#' @param save A logical value. If TRUE (default), saves the plot to file.
#' @return A ggplot object displaying the data trends.
#' @export
plot_data_trends_all_outcome <- function(model,
                             display_ci = TRUE,
                             recenter_t = TRUE,
                             display_vertical_event_lines = TRUE,
                             ylim_min = -0.1, ylim_max = 1.1,
                             save = TRUE) {
  # Extract outcome name and variable name
  empirical_spec <- model$empirical_spec
  outcome_name <- model$empirical_spec$outcome_name
  outcome_varname <- model$empirical_spec$outcome_varname
  x_lab <- model$empirical_spec$x_lab
  J <- length(model$priors)

  # Extract data
  y_g <- get_y_numeric_and_g_from_model(model)
  y <- y_g$y
  outcomes <- model$Y$names
  g <- y_g$g
  g_max <- max(g)

  plot_data_dfs <- lapply(outcomes,
    function(outcome) get_weighted_means_by_gjt(1 * (y == outcome), g) %>%
      mutate(outcome = outcome))
  plot_data_df <- do.call(rbind, plot_data_dfs)
  shapetypes <- get_shape_types(outcomes)
  shape_labels <- outcomes
  shape_manual <- ggplot2::scale_shape_manual(values = shapetypes,
                                              labels = shape_labels)

  p <- plot_data_df %>%
    dplyr::mutate(t = as.integer(t - ifelse(recenter_t, g_max - 1, 1))) %>%
    ggplot2::ggplot(ggplot2::aes_string(x = "t", y = "weighted_mean",
      group = "interaction(g, j, outcome)",
      color = "outcome",
      linetype = "g", linesize = "g", shape = "outcome")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    add_vertical_lines_at_g(unfactor_numerics(plot_data_df$g) - 1 -
      ifelse(recenter_t, (g_max - 1), 0),
      display = display_vertical_event_lines) +
    reference_hline() +
    ggplot2::scale_x_continuous(breaks = integer_breaks()) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_linetype_manual(values = get_line_types(g),
      labels = c("Control", "Treated")) +
    ggplot2::scale_color_manual(
      values = RColorBrewer::brewer.pal(max(3, length(outcomes)), "Set1")[1:length(outcomes)]) +
    ggplot2::labs(x = empirical_spec$x_lab,
      y = "Share of Individuals",
      color = "Outcome",
      linetype = "Treatment Status", shape = "Outcome")

  p_with_title <- p +
    ggplot2::ggtitle(empirical_spec_to_plot_title(
      empirical_spec, "Trends", display_se_info = FALSE
    ))

  # Save plots if requested
  if (save) {
    save_plot(get_plot_averages_all_outcome_directory(empirical_spec), p)
    save_plot(get_plot_averages_all_outcome_directory(empirical_spec,
      suffix = "-with-title"), p_with_title)
  }

  p_with_title
}



#' Plot Data Trends Conditional on Pretreatment Outcome
#'
#' Creates a plot showing outcome trends over time, conditional on having a specific pretreatment outcome value.
#'
#' @param model A model object containing outcome data, group information, treatment period, and empirical specifications.
#' @param conditional_y_value_pretreatment Numeric. The pretreatment outcome value to condition on. Default is 1.
#' @param display_ci Logical. Whether to display confidence intervals in the plot. Default is TRUE.
#' @param ylim_min Minimum value for the y-axis. Default is 0.
#' @param ylim_max Maximum value for the y-axis. Default is 1.
#' @return A ggplot object displaying the conditional data trends for the specified outcome variable.
#' @export
plot_data_trends_conditional <- function(model,
                                         conditional_y_value_pretreatment = 1,
                                         display_ci = TRUE,
                                         ylim_min = -0.1, ylim_max = 1.1) {
  # Extract outcome name and variable name
  empirical_spec <- model$empirical_spec
  outcome_name <- model$empirical_spec$outcome_name
  x_lab <- model$empirical_spec$x_lab

  # Construct data frame
  Y_names_of_interest <- model$Y$names[model$Y$values == conditional_y_value_pretreatment]
  y_g <- get_y_numeric_and_g_from_model(model)
  y_at_pretreatment <- generate_df_y_pretreatment(model) %>%
    dplyr::pull(y_at_pretreatment)
  y <- y_g$y_numeric[y_at_pretreatment %in% Y_names_of_interest,]
  g <- y_g$g[y_at_pretreatment %in% Y_names_of_interest]


  # Plot data
  weighted_means_by_gjt <- get_weighted_means_by_gjt(y, g)
  p <- plot_weighted_averages(weighted_means_by_gjt,
                              display_ci = display_ci) +
    ggplot2::scale_color_manual(values = "Black", name = "", guide="none") +
    ylim(0, 1) +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(outcome_name) +
    ggplot2::ylim(ylim_min, ylim_max)

  # Save plots
  p_with_title <- p +
    ggplot2::ggtitle(empirical_spec_to_plot_title(
      empirical_spec, paste0("Trends Conditional on Pretreatment y = ",
                             conditional_y_value_pretreatment),
      display_se_info = FALSE # Not cluster se
    ))
  save_plot(get_plot_conditional_averages_directory(empirical_spec,
    suffix = paste0("-", conditional_y_value_pretreatment)), p)
  save_plot(get_plot_conditional_averages_directory(empirical_spec,
    suffix = paste0("-", conditional_y_value_pretreatment, "-with-title")), p_with_title)

  p_with_title
}


#' Save Plot
#'
#' Saves a ggplot object to the specified file, creating the directory if needed.
#'
#' @param filename A character string specifying the file transition to save the plot.
#' @param plot_to_be_saved The ggplot object to be saved.
#' @param width Numeric. Width of the saved plot in inches.
#' Default is get_empirical_plot_default_width().
#' @param height Numeric. Height of the saved plot in inches.
#' Default is get_empirical_plot_default_height().
#' @return None. The plot is saved to the specified file.
#' @export
save_plot <- function(filename, plot_to_be_saved, width = NULL, height = NULL) {
  dir_name <- dirname(filename)
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
  }
  if (is.null(width)) {
    width <- get_empirical_plot_default_width()
  }
  if (is.null(height)) {
    height <- get_empirical_plot_default_height()
  }

  ggsave(filename, plot_to_be_saved, width = width, height = height)
}

#' Get Default Width for Empirical Plots
#'
#' Returns the default width (in inches) used for empirical plots.
#'
#' @return A numeric value indicating the default width. Default is 8.
#' @examples
#' get_empirical_plot_default_width()
#' @export
get_empirical_plot_default_width <- function() {
  return(8)
}

#' Get Default Height for Empirical Plots
#'
#' Returns the default height (in inches) used for empirical plots.
#'
#' @return A numeric value indicating the default height. Default is 5.
#' @examples
#' get_empirical_plot_default_height()
#' @export
get_empirical_plot_default_height <- function() {
  return(5)
}


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

#' Generate ATT Decomposition Plots for a Given Outcome
#'
#' Produces and saves two decomposition plots (unweighted and weighted)
#' for a specified outcome using the supplied causal model.
#' Returns the generated plots with informative titles.
#'
#' @param outcome_of_interest A string specifying the outcome variable for which
#' decomposition plots should be generated.
#' @param model A transition model estimate
#'
#' @return A list of two ggplot objects:
#' \itemize{
#'   \item{ATT Decomposition by Flow (unweighted)}
#'   \item{ATT Decomposition by Flow (weighted)}
#' }
#' The plots are also saved as PDF files in a subdirectory under `plots/`.
#'
#' @examples
#' \dontrun{
#'   model <- estimate_causal_model(...)  # Assume this creates a suitable model
#'   generate_decomposition_plot_pair("employment", model)
#' }
#'
#' @param outcome_label_function Optional function to format outcome labels in
#'   the legend. Takes a character vector and returns formatted labels.
#' @export
generate_decomposition_plot_pair <- function(outcome_of_interest, model,
                                             outcome_label_function = NULL) {
  empirical_spec <- model$empirical_spec
  empirical_spec$outcome_name <- outcome_of_interest # for plot title
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name
  x_lab <- model$empirical_spec$x_lab

  # Get number of outcomes for color palette
  n_outcomes <- length(model$Y$names)
  # If outcome_label_function is provided, use formatted names for palette
  palette_names <- if (!is.null(outcome_label_function)) {
    outcome_label_function(model$Y$names)
  } else {
    model$Y$names
  }
  color_palette <- RColorBrewer::brewer.pal(max(3, n_outcomes), "Set1")[1:n_outcomes]
  names(color_palette) <- palette_names

  p_catt <- plot_decomposition_catt(model, outcome_of_interest,
                                    outcome_label_function = outcome_label_function) +
    scale_color_manual(values = color_palette) +
    ggplot2::xlab(x_lab)
  p_catt_weighted <- plot_decomposition_catt_weighted(model, outcome_of_interest,
                                                      outcome_label_function = outcome_label_function) +
    scale_color_manual(values = color_palette) +
    ggplot2::xlab(x_lab)

  save_plot(paste0("plots/", outcome_varname, "/",
                   outcome_name, "-decomposition-by-catt.pdf"),
            p_catt)
  save_plot(paste0("plots/", outcome_varname, "/",
                   outcome_name, "-decomposition-by-catt-weighted.pdf"),
            p_catt_weighted)
  list(p_catt +
         ggplot2::ggtitle(empirical_spec_to_plot_title(empirical_spec,
                                                       "ATT Decomposition by Flow")),
       p_catt_weighted +
         ggplot2::ggtitle(empirical_spec_to_plot_title(empirical_spec,
                                                       "ATT Decomposition by Flow (Weighted)")))
}

#' Generate ATT Decomposition Plots for All Outcomes
#'
#' Iterates over all outcomes in the model and generates/saves decomposition plots
#' for each using \code{generate_decomposition_plot_pair}.
#'
#' @param model A fitted model object containing multiple outcome specifications.
#'
#' @return A list of lists of ggplot objects, each corresponding to a pair of
#' decomposition plots (unweighted and weighted) for each outcome.
#'
#' @examples
#' \dontrun{
#'   model <- estimate_causal_model(...)
#'   generate_all_decomposition_plots(model)
#' }
#'
#' @export
generate_all_decomposition_plots <- function(model, outcome_label_function = NULL) {
  lapply(model$Y$names, generate_decomposition_plot_pair, model,
         outcome_label_function = outcome_label_function)
}

#' Generate Per-Type Decomposition Plot Pair
#'
#' Generates and saves decomposition plots faceted by latent type, showing
#' type-specific weighted CATTs. Follows the pattern of
#' \code{generate_decomposition_plot_pair} but uses \code{plot_decomposition_by_type}.
#'
#' @param outcome_of_interest The outcome variable of interest.
#' @param model A fitted TransitionModel object.
#'
#' @return A list of two ggplot objects (unweighted and weighted by-type decomposition).
#'   Plots are saved as PDF files.
#'
#' @param outcome_label_function Optional function to format outcome labels in
#'   the legend. Takes a character vector and returns formatted labels.
#' @export
generate_decomposition_plot_pair_by_type <- function(outcome_of_interest, model,
                                                      outcome_label_function = NULL) {
  empirical_spec <- model$empirical_spec
  empirical_spec$outcome_name <- outcome_of_interest
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name
  x_lab <- model$empirical_spec$x_lab
  J <- length(model$priors)

  n_outcomes <- length(model$Y$names)
  palette_names <- if (!is.null(outcome_label_function)) {
    outcome_label_function(model$Y$names)
  } else {
    model$Y$names
  }
  color_palette <- RColorBrewer::brewer.pal(max(3, n_outcomes), "Set1")[1:n_outcomes]
  names(color_palette) <- palette_names

  p_catt <- plot_decomposition_by_type(model, outcome_of_interest, y_to_use = "catt",
                                       outcome_label_function = outcome_label_function) +
    scale_color_manual(values = color_palette) +
    ggplot2::xlab(x_lab)
  p_catt_weighted <- plot_decomposition_by_type(model, outcome_of_interest, y_to_use = "catt_weighted",
                                                outcome_label_function = outcome_label_function) +
    scale_color_manual(values = color_palette) +
    ggplot2::xlab(x_lab)

  save_plot(paste0("plots/", outcome_varname, "/",
                   outcome_name, "-decomposition-by-type-J", J, "-catt.pdf"),
            p_catt)
  save_plot(paste0("plots/", outcome_varname, "/",
                   outcome_name, "-decomposition-by-type-J", J, "-catt-weighted.pdf"),
            p_catt_weighted)
  list(p_catt +
         ggplot2::ggtitle(empirical_spec_to_plot_title(empirical_spec,
                                                       paste0("ATT Decomposition by Flow (J=", J, ")"))),
       p_catt_weighted +
         ggplot2::ggtitle(empirical_spec_to_plot_title(empirical_spec,
                                                       paste0("ATT Decomposition by Flow, Weighted (J=", J, ")"))))
}

#' Generate Per-Type Decomposition Plots for All Outcomes
#'
#' Iterates over all outcomes in the model and generates/saves per-type
#' decomposition plots using \code{generate_decomposition_plot_pair_by_type}.
#'
#' @param model A fitted TransitionModel object with J > 1.
#' @return A list of lists of ggplot objects.
#' @export
generate_all_decomposition_plots_by_type <- function(model,
                                                      outcome_label_function = NULL) {
  lapply(model$Y$names, generate_decomposition_plot_pair_by_type, model,
         outcome_label_function = outcome_label_function)
}

