#' Plot Expected Outcomes from Transition Model
#'
#' Generates a plot of the expected outcomes for control and treated groups across all time periods
#' as defined in the transition model. This visualization aids in understanding the dynamic effects
#' of treatment over time and across different groups.
#'
#' @param model A list containing the transition model elements, including priors, transition matrices
#'        (`Ps_treated` and `Ps_control`), and outcome values (`Y$values`). The transition model should also
#'        include `treated_period`, indicating the period in which treatment begins, and `T_max`,
#'        the maximum time period considered in the model.
#' @param y_name A scalar value of y to compute transition probabilities
#' @param y_past_name A scalar value of y in the past to compute transition probabilities
#' @param aggregate Logical; if `TRUE`, the plot aggregates outcomes over all groups (j), showing a summary
#'        outcome for each time period (`t`) and treatment status (`g`). If `FALSE`, outcomes are plotted
#'        for each group and time period separately.
#' @param display_vertical_event_lines whether to display vertical lines for events
#' @param recenter_t whether to recenter the time axis (g is displayed as 0)
#' @param y_lab the label for the y-axis
#' @param display_ci Logical; if `TRUE` (default), display confidence intervals.
#' @param uniform_ci Logical; if `TRUE` (default), compute uniform (simultaneous)
#'        confidence intervals; if `FALSE`, compute pointwise intervals.
#' @return a ggplot object
#' @export
plot_transition_model_transition_probabilities <- function(model, y_name, y_past_name, aggregate = FALSE,
                                                           display_vertical_event_lines = TRUE,
                                                           recenter_t = TRUE,
                                                           y_lab = NULL,
                                                           display_ci = TRUE,
                                                           uniform_ci = TRUE) {
  # Get indices
  current_index <- which(model$Y$names == y_name)[1]
  past_index <- which(model$Y$names == y_past_name)[1]

  # Return errors if the indices are not found
  if (is.na(current_index)) {
    stop("y_name not found in model outcome space (model$Y)")
  }
  if (is.na(past_index)) {
    stop("y_past_name not found in model outcome space (model$Y)")
  }

  # Get expected outcomes
  plot_data_df <- get_transition_probability_df(model, aggregate = aggregate) %>%
    dplyr::filter(y_current == y_name, y_past == y_past_name) %>%
    dplyr::rename(y = probability) %>%
    attach_ci_to_plot_data_df(model,
                              "Ps_control_empirical",
                              "Ps_treated_empirical",
                              current_index = current_index,
                              past_index = past_index,
                              uniform_ci = uniform_ci) %>%
    dplyr::mutate(g = as.factor(g), j = factor_j_with_priors(j, model$priors))



  # Plot them
  if (is.null(y_lab)) {
    y_value_name <- ifelse(y_name == 0, "Unemployment",
                           ifelse(y_name == 1, "Employment",
                                  ifelse(y_name == -1, "OLF", y_name)))
    y_value_past_name <- ifelse(y_past_name == 0, "unemployed",
                                ifelse(y_past_name == 1, "employed",
                                       ifelse(y_past_name == -1, "OLF", y_past_name)))
    y_lab <- paste0(y_value_name,
                    " probability conditional on past (",
                    y_value_past_name, ")")
  }


  plot_from_plot_data_df(plot_data_df,
                         recenter_t = recenter_t,
                         display_vertical_event_lines = display_vertical_event_lines,
                         y_lab = y_lab)

}


#' Plot Transition Probability Differences from Transition Model
#'
#' Generates a plot of the differences in transition probabilities between
#' treated and control groups across all time periods as defined in the
#' transition model. This visualization aids in understanding the treatment
#' effects on transition probabilities over time.
#'
#' @param model A list containing the transition model elements, including priors, transition matrices
#'        (`Ps_treated` and `Ps_control`), and outcome values (`Y$values`). The transition model should also
#'        include `treated_period`, indicating the period in which treatment begins, and `T_max`,
#'        the maximum time period considered in the model.
#' @param y_name A scalar value of y to compute transition probabilities
#' @param y_past_name A scalar value of y in the past to compute transition probabilities
#' @param aggregate Logical; if `TRUE`, the plot aggregates outcomes over all groups (j), showing a summary
#'        outcome for each time period (`t`) and treatment status (`g`). If `FALSE`, outcomes are plotted
#'        for each group and time period separately.
#' @param display_vertical_event_lines whether to display vertical lines for events
#' @param recenter_t whether to recenter the time axis (g is displayed as 0)
#' @param y_lab the label for the y-axis
#' @param display_ci Logical; if `TRUE` (default), display confidence intervals.
#' @param uniform_ci Logical; if `TRUE` (default), compute uniform (simultaneous)
#'        confidence intervals; if `FALSE`, compute pointwise intervals.
#' @return a ggplot object
#' @export
plot_transition_model_transition_probability_differences <- function(model, y_name, y_past_name, aggregate = FALSE,
                                                           display_vertical_event_lines = TRUE,
                                                           recenter_t = TRUE,
                                                           y_lab = NULL,
                                                           display_ci = TRUE,
                                                           uniform_ci = TRUE) {
  # Get indices
  current_index <- which(model$Y$names == y_name)[1]
  past_index <- which(model$Y$names == y_past_name)[1]

  # Return errors if the indices are not found
  if (is.na(current_index)) {
    stop("y_name not found in model outcome space (model$Y)")
  }
  if (is.na(past_index)) {
    stop("y_past_name not found in model outcome space (model$Y)")
  }

  # Get expected outcomes
  treated_period <- model$treated_period
  J <- length(model$priors)
  ts <- 2:(length(model$Ps_control[[1]]) + 1)
  summary <- summarize_transition_model(model)
  plot_data_df_list <- lapply(1:J, function(j) {
    y <- fetch_parameter(summary, "Ps_empirical_diffs", j, current_index, past_index)
    data.frame(t = ts, y = y, g = treated_period, j = j) }
  )
  plot_data_df <- do.call(rbind, plot_data_df_list) %>%
    attach_ci_to_plot_data_df(model,
                              "Ps_empirical_diffs",
                              current_index = current_index,
                              past_index = past_index,
                              uniform_ci = uniform_ci)

  # Plot them
  if (is.null(y_lab)) {
    y_value_name <- ifelse(y_name == 0, "Unemployment",
                           ifelse(y_name == 1, "Employment",
                                  ifelse(y_name == -1, "OLF", y_name)))
    y_value_past_name <- ifelse(y_past_name == 0, "unemployed",
                                ifelse(y_past_name == 1, "employed",
                                       ifelse(y_past_name == -1, "OLF", y_past_name)))
    y_lab <- paste0(y_value_name,
                    " probability conditional on past (",
                    y_value_past_name, ")")
  }

  plot_from_plot_data_df(plot_data_df,
                         recenter_t = recenter_t,
                         display_vertical_event_lines = display_vertical_event_lines,
                         y_lab = y_lab) +
    reference_hline()

}

#' Plot treatment effects by transition pattern group
#' @export
#' @param model A TransitionModel object to be plotted
#' @param display_ci Whether to display confidence intervals if available (default is FALSE)
#' @param display_vertical_event_lines whether to display vertical lines for events
#' @param recenter_t whether to recenter the time axis (g is displayed as 0)
#' @param uniform_ci Logical; if `TRUE` (default), compute uniform (simultaneous)
#'        confidence intervals; if `FALSE`, compute pointwise intervals.
#' @return a ggplot object
plot_ltatt <- function(model, display_ci = FALSE,
                        display_vertical_event_lines = TRUE,
                        recenter_t = TRUE,
                        uniform_ci = TRUE) {
  # Data frame for LTATTs
  J <- length(model$priors)
  treated_period <- model$treated_period
  pretreatment_periods <- 1:(treated_period - 1)
  T_max <- length(model$Ps_control[[1]]) + 1
  model_summary <- summarize_transition_model(model)
  ltatt <- sapply(model_summary$LTATTs, c)
  gatts_matrix <- rbind(matrix(NA, nrow = length(pretreatment_periods), ncol = J),
                        ltatt)

  # Grab gs and js
  gs <- (c(max(pretreatment_periods) + 1))
  js <- (1:J)
  ts <- (1:T_max)

  # Make into long data frame
  plot_data_df <- data.frame(y = as.vector(gatts_matrix),
                         j = rep(js, each = length(ts)),
                         g = (rep(gs[1], length(ts)*J)),
                         t = rep(ts, J)) %>%
    filter(t >= treated_period) %>%
    attach_ci_to_plot_data_df(model, "LTATTs", uniform_ci = uniform_ci) %>%
    dplyr::mutate(j = factor_j_with_priors(j, model$priors))

  # Attach confidence intervals
  p <- plot_from_plot_data_df(plot_data_df,
                         recenter_t = recenter_t,
                         display_vertical_event_lines = display_vertical_event_lines,
                         display_lines = !display_ci)

  p +
    reference_hline() +
    ggplot2::ylab("Latent Type Average Treatment Effects on Treated (LTATT)")
}

#' Generate Plot Data for LTATTs by Pretreatment Outcomes
#'
#' Computes the Latent Type Average Treatment Effects on the Treated (LTATTs) by pretreatment outcomes.
#' Within each `initial_outcome`, `j`, and `t`, the function calculates the difference
#' in `y` by `g`.
#'
#' @param model A fitted model object used for counterfactual analysis.
#' @param ts A time series or panel data structure required for generating counterfactuals.
#'
#' @return A data frame containing the computed LTATTs, with columns:
#'   - `initial_outcome`: The initial outcome variable.
#'   - `j`: The grouping variable.
#'   - `t`: The time period.
#'   - `y`: The computed difference in outcomes.
#'   - `type`: A combined factor of `j` and `initial_outcome`.
#'
#' @details This function first generates counterfactual data using
#'   `generate_plot_data_df_for_counterfactuals()`, then groups the data by
#'   `initial_outcome`, `j`, and `t`, computes the difference in `y`, and
#'   finally adds a `type` column by interacting `j` and `initial_outcome`.
#'
#' @seealso \code{\link{generate_plot_data_df_for_counterfactuals}}
#'
#' @export
generate_plot_data_df_for_ltatt_by_pretreatment_outcomes <- function(model, ts) {
  # Compute LTATTs by pretreatment outcomes; within each initial_outcome, j, t,
  # Take difference of y by g
  generate_plot_data_df_for_counterfactuals(model, ts) %>%
    group_by(initial_outcome, j, t) %>%
    summarise(y = diff(y)) %>%
    mutate(type = interaction(j, initial_outcome))
}

#' Plot LTATTs by pretreatment outcomes
#'
#' This function plots latent group average treatment effects on treated (LTATTs)
#' by pretreatment outcomes, displaying treatment effect heterogeneity across
#' different initial outcomes just before treatment.
#'
#' @export
#' @param model A TransitionModel object to be plotted
#' @param weighted A logical value. If TRUE, display LTATT * pre-treatment pmfs. Default is FALSE.
#' @return A ggplot object visualizing LTATTs by pretreatment outcomes
plot_ltatt_by_pretreatment_outcomes <- function(model, weighted = FALSE) {
  # Extract data
  J <- length(model$priors)
  treated_period <- model$treated_period
  # Use Ps_control_empirical for consistency with generate_plot_data_df functions
  ts_post_treatment <- (treated_period-1):(length(model$Ps_control_empirical[[1]]))
  ts <- c(treated_period-1, (ts_post_treatment+1))
  pmf_pre_treatment_treated <- model$pmfs_pre_treatment_treated[[1]]

  # Generate plot data for counterfactuals
  plot_data_df <- generate_plot_data_df_for_ltatt_by_pretreatment_outcomes(model, ts) %>%
    mutate(linetype = initial_outcome) %>%
    mutate(g = treated_period) %>%
    # Extract outcome name from labeled factor for matching
    mutate(weights = pmf_pre_treatment_treated[
      match(extract_outcome_name(initial_outcome), model$Y$names)
    ]) %>%
    mutate(catt = y, catt_weighted = catt * weights) %>%
    mutate(j = factor_j_with_priors(j, model$priors))

  # Plot counterfactual outcomes, by initial outcome and `j` pairs using ggplot2
  # Remove warnings for scale already present message due to scales attached in plot_from_plot_data_df
  shapetypes <- get_shape_types_for_unique_data(unique(plot_data_df$linetype))
  suppressMessages(
    p <- plot_from_plot_data_df(plot_data_df, point_size = 3) +
      reference_hline() +
      ggplot2::labs(x = "Time", y = "Latent Type Average Treatment Effects on Treated (LTATT)",
                    color = "Latent Type",
                    shape = "Pretreatment Outcome") +
      ggplot2::scale_shape_manual(values=shapetypes, labels=unique(plot_data_df$linetype)) +
      ggplot2::scale_linetype_manual(values=rep("solid", length(shapetypes))) +
      ggplot2::guides(linetype = "none",
                      color = if (J == 1) "none" else ggplot2::guide_legend())
  )
  p
}

#' Plot Event Study Results from a Data Frame
#'
#' This function generates a plot from a data frame containing event study results.
#' The plot includes points, lines, and optional confidence intervals. It supports
#' customization of axis labels, vertical event lines, and recentering of time periods.
#'
#' @param plot_data_df A data frame containing the event study data. Must include columns `t` (time),
#' `y` (outcome), `g` (group), and `j` (cohort). Optionally, it can include columns `ci_lower` and `ci_upper` for confidence intervals.
#' @param recenter_t Logical. If `TRUE`, recenters the time variable `t` relative to the maximum value of `g`.
#' @param display_vertical_event_lines Logical. If `TRUE`, adds vertical dashed lines for event time points.
#' @param display_lines Logical. If `TRUE`, displays lines connecting points.
#' @param y_lab A character string specifying the label for the y-axis. Defaults to `NULL`, in which case a default label is applied.
#' @param point_size A numeric value specifying the size of the points. Defaults to `1`.
#'
#' @return A ggplot object representing the event study plot.
#'
#' @details
#' - The function uses `ggplot2` for visualization and assumes the input data frame has been preprocessed for plotting.
#' - If the `type` column is not already present in `plot_data_df`, it is created as the interaction of `g` and `j`.
#' - Confidence intervals are displayed using `geom_errorbar` if the `display_ci` argument is `TRUE` and the relevant columns are available.
#' - Vertical event lines are controlled using the `display_vertical_event_lines` argument and added at appropriate time points.
#'
#' @examples
#' # Example data frame
#' plot_data_df <- data.frame(
#'   t = c(-1, 0, 1, -1, 0, 1),
#'   y = c(0.2, 0.3, 0.4, 0.1, 0.15, 0.2),
#'   g = c(0, 0, 0, 1, 1, 1),
#'   j = c(1, 1, 1, 2, 2, 2),
#'   ci_lower = c(0.1, 0.2, 0.3, 0.05, 0.1, 0.15),
#'   ci_upper = c(0.3, 0.4, 0.5, 0.15, 0.2, 0.25)
#' )
#'
#' # Generate the plot
#' plot_from_plot_data_df(
#'   plot_data_df,
#'   display_ci = TRUE,
#'   recenter_t = TRUE,
#'   display_vertical_event_lines = TRUE,
#'   y_lab = "Outcome"
#' )
#'
#' @import ggplot2
#' @importFrom dplyr mutate
plot_from_plot_data_df <- function(plot_data_df,
                                   recenter_t = TRUE,
                                   display_vertical_event_lines = TRUE,
                                   display_lines = TRUE,
                                   y_lab = NULL,
                                   point_size = NULL) {
  # Plot them
  g_max <- max(unfactor_numerics(plot_data_df$g))

  # Check if there is a separate column for type; if not, use interaction(g,j)
  if (!("type" %in% names(plot_data_df))) {
    plot_data_df <- plot_data_df %>%
      dplyr::mutate(type = interaction(g, j))
  }
  if (!("linetype" %in% names(plot_data_df))) {
    plot_data_df <- plot_data_df %>%
      dplyr::mutate(linetype = g)
  }
  display_ci <- FALSE
  if (("ci_lower" %in% names(plot_data_df)) && ("ci_upper" %in% names(plot_data_df))) {
    display_ci <- TRUE
  }
  # Create a single position_dodge object to ensure consistent alignment
  dodge_width <- ifelse(display_ci, 0.5, 0)
  pd <- ggplot2::position_dodge(width = dodge_width)
  shape_varname <- ifelse("initial_outcome" %in% names(plot_data_df),
                          "initial_outcome", "linetype")

  # Convert to factors, preserving existing factor levels/labels if already set
  p <- plot_data_df %>%
    dplyr::mutate(
      g = if (is.factor(g)) g else as.factor(g),
      j = if (is.factor(j)) j else as.factor(j),
      linetype = if (is.factor(linetype)) linetype else as.factor(linetype)
    ) %>%
    dplyr::mutate(t = as.integer(t-ifelse(recenter_t, g_max-1, 1))) %>% # To account for t = 0
    ggplot2::ggplot(ggplot2::aes_string(x = "t", y = "y",
                                        group = "type",
                                        color = "j", linetype = "linetype", linesize = "g",
                                        shape = shape_varname)) +
    ggplot2::geom_point(size = ifelse(is.null(point_size), 2, point_size),
                        position = pd)


  # Add confidence intervals if required
  if (display_ci) {
    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                             width = 0.5, linetype = "solid",
                             position = pd)
    display_lines <- FALSE
  }

  if (display_lines) {
    p <- p +
      ggplot2::geom_line()
  }

  # Add theme and vertical event lines
  p <- p +
    get_transition_plot_theme(plot_data_df$linetype,
                              plot_data_df$j,
                              plot_data_df[,p$labels$shape],
                              y_lab) +
    add_vertical_lines_at_g(unfactor_numerics(plot_data_df$g) - 1 - ifelse(recenter_t, (g_max-1), 0),
                            display = display_vertical_event_lines) +
    ggplot2::ylab(y_lab) +
    ggplot2::scale_x_continuous(breaks = integer_breaks())
  p
}

#' Plot Counterfactual Comparisons by Pretreatment Outcomes
#'
#' This function generates plots showing the expected outcomes for treated and untreated groups,
#' based on pretreatment outcome values and transition models.
#'
#' @param model A list representing the transition model.
#' @return A \code{ggplot} object displaying expected outcomes for treated and untreated groups,
#'         stratified by initial outcome values and transition model indices.
#'
#' @details
#' The function calculates expected outcomes for control (untreated) and treated scenarios,
#' starting from different initial outcome values. It iterates over all initial values and transition
#' model groups, then combines results into a single data frame for plotting.
#'
#' @examples
#' # Example usage:
#' model <- generate_transition_model()
#' plot_counterfactuals_by_pretreatment_outcomes(model)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
plot_counterfactuals_by_pretreatment_outcomes <- function(model,
                                                          diff_by_status = TRUE,
                                                          uniform_ci = TRUE) {
  # Extract data
  J <- length(model$priors)
  treated_period <- model$treated_period
  # Use Ps_control_empirical for ts since that's what generate_plot_data_df_for_counterfactuals uses
  ts <- 1:(length(model$Ps_control_empirical[[1]]) + 1)

  # Generate plot data
  plot_data_df <- generate_plot_data_df_for_counterfactuals(model, ts,
                                                            diff_by_status = diff_by_status)

  if (!is.null(model$bootstrap_estimates)) {
    estimate <- plot_data_df$y
    bootstrap_estimates <- lapply(model$bootstrap_estimates,
                                  generate_plot_data_df_for_counterfactuals, ts,
                                  diff_by_status = diff_by_status) %>% lapply(pull, y)
    plot_data_df <- bootstrap_estimates_to_ci_df(estimate, bootstrap_estimates, uniform_ci = uniform_ci) %>%
      select(-t) %>%
      cbind(plot_data_df) %>%
      arrange(j)

  }

  # Plot counterfactual outcomes, by initial outcome and `j` pairs using ggplot2
  p <- plot_from_plot_data_df(plot_data_df) +
    ggplot2::labs(shape = "Pretreatment Outcome")

  save_plot(get_plot_counterfactuals_by_pretreatment_outcomes_directory(model$empirical_spec, J),
            p)

  p +
    ggplot2::ggtitle(empirical_spec_to_plot_title(
      model$empirical_spec, "Counterfactuals by Pretreatment Outcomes"
    ))
}

#' Plot Counterfactual Differences by Pretreatment Outcomes
#'
#' Plots the difference between treated and control counterfactual outcomes
#' (treated - control) for each pretreatment outcome and latent group,
#' with confidence interval bands.
#'
#' @param model A fitted transition model with bootstrap_estimates for CIs.
#' @param uniform_ci Logical; if TRUE, use uniform CIs. Default TRUE.
#' @param recenter_t Logical; if TRUE, recenter time so treatment period is 0.
#'
#' @return A ggplot object showing differences with CI bands.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
plot_counterfactuals_differences_by_pretreatment_outcomes <- function(
    model,
    uniform_ci = TRUE,
    recenter_t = TRUE) {

  # Extract data
  J <- length(model$priors)
  treated_period <- model$treated_period
  # Use Ps_control_empirical for consistency with generate_plot_data_df_for_counterfactuals_differences
  ts_post_treatment <- (treated_period - 1):(length(model$Ps_control_empirical[[1]]))
  ts <- c(treated_period - 1, (ts_post_treatment + 1))

 # Generate differences data with CIs
  plot_data_df <- generate_plot_data_df_for_counterfactuals_differences(
    model, ts, uniform_ci = uniform_ci
  )

  # Recenter time if requested
  # 0 = last pre-treatment period, 1 = first post-treatment period
  if (recenter_t) {
    plot_data_df <- plot_data_df %>%
      mutate(t = t - (treated_period - 1))
  }

  # Filter out t=0 (the conditioning period) where difference is always 0 by construction
  # since we condition on the pretreatment outcome
  plot_data_df <- plot_data_df %>%
    filter(t >= ifelse(recenter_t, 1, treated_period))

  # Check if CIs are available
  has_ci <- !all(is.na(plot_data_df$ci_lower))

  # Create interaction for dodging (combination of initial_outcome and j)
  plot_data_df <- plot_data_df %>%
    mutate(dodge_group = interaction(initial_outcome, j))

  # Get colors from theme
  linecolors <- get_line_colors(unique(plot_data_df$j))

  # Position dodge for spacing between points at same time
  pd <- ggplot2::position_dodge(width = 0.3)

  # Build plot - color by j, shape by initial_outcome
  p <- ggplot2::ggplot(plot_data_df, ggplot2::aes(
    x = t,
    y = estimate,
    color = j,
    shape = initial_outcome,
    group = dodge_group
  )) +
    ggplot2::geom_point(size = 3, position = pd)

  # Add CI error bars if available
  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2,
      position = pd
    )
  }

  # Add reference lines and styling
  p <- p +
    reference_hline() +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(values = linecolors) +
    ggplot2::labs(
      x = "Year Relative to the Policy Change",
      y = "Difference (Treated - Control)",
      color = "Latent Type",
      shape = "Pretreatment Outcome"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 1),
      color = if (J == 1) "none" else ggplot2::guide_legend(order = 2)
    )

  # Save plot
  save_plot(
    get_plot_counterfactuals_by_pretreatment_outcomes_directory(
      model$empirical_spec, J, suffix = "-differences"
    ),
    p
  )

  # Return with title
  p + ggplot2::ggtitle(
    empirical_spec_to_plot_title(
      model$empirical_spec,
      "Counterfactual Differences by Pretreatment Outcomes"
    )
  )
}

#' Generate Counterfactual Data for Given Outcomes and Groups
#'
#' @param Y_names Names of the outcomes.
#' @param ts Time periods to be used.
#' @param diff_by_status Logical; if `TRUE`, compute differences by treatment status.
#'
#' @return A combined data frame of counterfactual outcomes for all subgroups and initial outcomes.
generate_plot_data_df_for_counterfactuals <- function(model, ts, diff_by_status = FALSE) {
  # Initialize an empty list to store rows
  treated_period <- model$treated_period
  Y_names <- model$Y$names
  Y_values <- model$Y$values
  J <- length(model$priors)
  ts_for_Ps <- ts[2:length(ts)]

  # Check if the length of time periods matches the length of transition probabilities
  if (length(model$Ps_control) != length(model$Ps_treated)) {
    stop("Length of control and treated transition probabilities must match.")
  }

  # Fetch Ps_control and Ps_treated for periods of interest
  Ps_control <- lapply(1:J, function(j) model$Ps_control_empirical[[j]])
  Ps_treated <- lapply(1:J, function(j) model$Ps_treated_empirical[[j]])

  # Iterate over initial outcomes and subgroups
  rows <- list()
  for (y_initial_index in 1:length(Y_names)) {
    for (j in 1:J) {
      # Set initial outcome vector
      initial_value <- numeric(length(Y_names))
      initial_value[y_initial_index] <- 1
      initial_outcome_name <- Y_names[y_initial_index]

      # Compute expected outcomes under control and treated scenarios
      expected_outcomes_control <- calculate_expected_outcomes_recursively(
        c(initial_value), Ps_control[[j]], Y_values
      )
      expected_outcomes_treated <- calculate_expected_outcomes_recursively(
        c(initial_value), Ps_treated[[j]], Y_values
      )

      # Create full y vectors (initial value + expected outcomes at each time)
      # y[t] gives the expected outcome at time t
      full_y_control <- c(Y_values[y_initial_index], expected_outcomes_control)
      full_y_treated <- c(Y_values[y_initial_index], expected_outcomes_treated)

      # Append rows for control and treated outcomes
      # Subset y to match the time periods in ts
      df_control <- data.frame(
        initial_outcome = initial_outcome_name, g = 0, j = j, t = ts,
        y = full_y_control[ts]
      )
      df_treated <- data.frame(
        initial_outcome = initial_outcome_name, g = treated_period, j = j, t = ts,
        y = full_y_treated[ts]
      )

      if (diff_by_status) {
        df_diff <- df_treated
        df_diff$y <- df_treated$y - df_control$y
        rows <- append(rows, list(df_diff))
      } else {
        rows <- append(rows, list(df_control)) %>% append(list(df_treated))
      }
    }
  }

  # Compute outcome shares for treated units at pretreatment period
  outcome_shares <- compute_outcome_shares(model, pmf_type = "pretreatment")

  # Combine all rows into a single data frame
  plot_data_df <- do.call(rbind, rows) %>%
    dplyr::mutate(j = factor_j_with_priors(j, model$priors),
                  initial_outcome = factor_outcome_with_share(
                    initial_outcome, Y_names, outcome_shares
                  )) %>%
    dplyr::mutate(type = paste(g, j, initial_outcome, sep = "_"))

  return(plot_data_df)
}

#' Generate Counterfactual Differences Data Frame (Treated - Control)
#'
#' Computes the difference between treated and control counterfactual outcomes
#' for each initial outcome and latent group, with optional bootstrap confidence
#' intervals (including uniform CIs).
#'
#' @param model A fitted transition model containing Ps_control_empirical,
#'   Ps_treated_empirical, and optionally bootstrap_estimates.
#' @param ts Time periods to compute counterfactuals for.
#' @param uniform_ci Logical; if TRUE, compute uniform (simultaneous) confidence
#'   intervals; if FALSE, compute pointwise intervals. Default is TRUE.
#'
#' @return A data frame with columns: initial_outcome, j, t, estimate, ci_lower,
#'   ci_upper (if bootstrap available), and type.
#'
#' @details
#' This function computes E(Y_t | Y_0, D=1) - E(Y_t | Y_0, D=0) for each initial
#' outcome Y_0 and latent group j. When bootstrap estimates are available,
#' confidence intervals are computed using the bootstrap distribution.
#'
#' @export
generate_plot_data_df_for_counterfactuals_differences <- function(model, ts,
                                                                   uniform_ci = TRUE) {
  # Extract model components
  treated_period <- model$treated_period
  Y_names <- model$Y$names
  Y_values <- model$Y$values
  J <- length(model$priors)
  bootstrap_models <- model$bootstrap_estimates

  # Fetch Ps for periods of interest

  ts_for_Ps <- ts[1:(length(ts) - 1)]

  # Helper function to compute differences for a given model

  compute_differences <- function(m) {
    Ps_control <- lapply(seq_len(J), function(j) m$Ps_control_empirical[[j]][ts_for_Ps])
    Ps_treated <- lapply(seq_len(J), function(j) m$Ps_treated_empirical[[j]][ts_for_Ps])

    # Use non-Markovian if available
    use_cumulative <- FALSE
    if (!is.null(m$Ps_control_from_pre) && !is.null(m$Ps_treated_from_pre)) {
      Ps_control <- lapply(seq_len(J), function(j) m$Ps_control_from_pre[[j]][ts_for_Ps])
      Ps_treated <- lapply(seq_len(J), function(j) m$Ps_treated_from_pre[[j]][ts_for_Ps])
      use_cumulative <- TRUE
    }

    diffs <- list()
    for (y_initial_index in seq_along(Y_names)) {
      for (j in seq_len(J)) {
        initial_value <- numeric(length(Y_names))
        initial_value[y_initial_index] <- 1

        # Use calculate_expected_outcomes for cumulative (non-Markovian) transitions,
        # or calculate_expected_outcomes_recursively for one-step Markov transitions
        if (use_cumulative) {
          expected_control <- calculate_expected_outcomes(
            initial_value, Ps_control[[j]], Y_values
          )
          expected_treated <- calculate_expected_outcomes(
            initial_value, Ps_treated[[j]], Y_values
          )
        } else {
          expected_control <- calculate_expected_outcomes_recursively(
            initial_value, Ps_control[[j]], Y_values
          )
          expected_treated <- calculate_expected_outcomes_recursively(
            initial_value, Ps_treated[[j]], Y_values
          )
        }

        # Difference: treated - control (starting from t=2 since t=1 is initial)
        diff_values <- expected_treated - expected_control
        key <- paste(Y_names[y_initial_index], j, sep = "_")
        diffs[[key]] <- diff_values
      }
    }
    diffs
  }

  # Compute point estimates
  point_estimates <- compute_differences(model)

  # Build result data frame
  rows <- list()
  for (y_initial_index in seq_along(Y_names)) {
    for (j in seq_len(J)) {
      key <- paste(Y_names[y_initial_index], j, sep = "_")
      estimate <- point_estimates[[key]]

      df_row <- data.frame(
        initial_outcome = Y_names[y_initial_index],
        j = j,
        t = ts[-1],  # Exclude first time point (initial outcome)
        estimate = estimate
      )

      # Add bootstrap CIs if available
      if (!is.null(bootstrap_models)) {
        bootstrap_diffs <- lapply(bootstrap_models, compute_differences)
        bootstrap_estimates <- lapply(bootstrap_diffs, function(d) d[[key]])

        ci_df <- bootstrap_estimates_to_ci_df(estimate, bootstrap_estimates,
                                              uniform_ci = uniform_ci)
        df_row$ci_lower <- ci_df$ci_lower
        df_row$ci_upper <- ci_df$ci_upper
      } else {
        df_row$ci_lower <- NA
        df_row$ci_upper <- NA
      }

      rows <- append(rows, list(df_row))
    }
  }

  # Compute outcome shares for treated units at pretreatment period
  outcome_shares <- compute_outcome_shares(model, pmf_type = "pretreatment")

  # Combine and format
  plot_data_df <- do.call(rbind, rows) %>%
    dplyr::mutate(
      j = factor_j_with_priors(j, model$priors),
      initial_outcome = factor_outcome_with_share(initial_outcome, Y_names, outcome_shares),
      type = paste(initial_outcome, j, sep = "_")
    )

  return(plot_data_df)
}

#' Generate Counterfactual Differences Data Using Transitions from Initial Period
#'
#' This function computes the difference (treated - control) in expected outcomes
#' using non-Markovian transition matrices from the initial period (t=0) to each
#' subsequent period. This is useful for analyzing treatment effects conditional
#' on the initial outcome state.
#'
#' @param model A fitted transition model object.
#' @param uniform_ci A logical; if TRUE, compute uniform confidence intervals.
#' @param include_post_treatment_periods A logical; if FALSE (default), only include
#'   pre-treatment periods and compute uniform CIs only on those periods.
#'
#' @return A data frame with columns: initial_outcome, j, t, estimate, ci_lower, ci_upper.
#' @export
generate_plot_data_df_for_counterfactuals_differences_initial_outcomes <- function(
    model, uniform_ci = TRUE, include_post_treatment_periods = FALSE) {

  # Extract model components
  Y_names <- model$Y$names
  Y_values <- model$Y$values
  J <- length(model$priors)
  bootstrap_models <- model$bootstrap_estimates
  treated_period <- model$treated_period

  # Check that Ps_from_0 are available
  if (is.null(model$Ps_control_from_0) || is.null(model$Ps_treated_from_0)) {
    stop("Model must have Ps_control_from_0 and Ps_treated_from_0. ",
         "Ensure estimate_post_em_params was called.")
  }

  # Time periods: from initial period (t=0) to the end
  T_max <- length(model$Ps_control_from_0[[1]]) + 1  # +1 because Ps has T-1 transitions

  # Determine which periods to include for uniform CI computation
  # t=1 is the initial/conditioning period, so pretreatment periods are t=2 to t=treated_period-1
  if (include_post_treatment_periods) {
    periods_for_ci <- 2:(T_max - 1)  # All periods after initial (excluding t=1)
  } else {
    periods_for_ci <- 2:(treated_period - 1)  # Only pre-treatment periods (excluding initial t=1)
  }

  # Helper function to compute differences for a given model
  compute_differences <- function(m) {
    Ps_control_from_0 <- m$Ps_control_from_0
    Ps_treated_from_0 <- m$Ps_treated_from_0

    diffs <- list()
    for (y_initial_index in seq_along(Y_names)) {
      for (j in seq_len(J)) {
        initial_value <- numeric(length(Y_names))
        initial_value[y_initial_index] <- 1

        # Use calculate_expected_outcomes (not recursively) because Ps_from_0
        # are already t-step cumulative transitions from the initial period,
        # not one-step Markov transitions
        expected_control <- calculate_expected_outcomes(
          initial_value, Ps_control_from_0[[j]], Y_values
        )
        expected_treated <- calculate_expected_outcomes(
          initial_value, Ps_treated_from_0[[j]], Y_values
        )

        # Difference: treated - control
        diff_values <- expected_treated - expected_control
        key <- paste(Y_names[y_initial_index], j, sep = "_")
        diffs[[key]] <- diff_values
      }
    }
    diffs
  }

  # Compute point estimates
  point_estimates <- compute_differences(model)

  # Build result data frame
  rows <- list()
  for (y_initial_index in seq_along(Y_names)) {
    for (j in seq_len(J)) {
      key <- paste(Y_names[y_initial_index], j, sep = "_")
      estimate_full <- point_estimates[[key]]

      # Filter to relevant periods
      estimate <- estimate_full[periods_for_ci]

      df_row <- data.frame(
        initial_outcome = Y_names[y_initial_index],
        j = j,
        t = periods_for_ci,
        estimate = estimate
      )

      # Add bootstrap CIs if available
      if (!is.null(bootstrap_models)) {
        bootstrap_diffs <- lapply(bootstrap_models, compute_differences)
        # Filter bootstrap estimates to relevant periods
        bootstrap_estimates <- lapply(bootstrap_diffs, function(d) d[[key]][periods_for_ci])

        ci_df <- bootstrap_estimates_to_ci_df(estimate, bootstrap_estimates,
                                              uniform_ci = uniform_ci)
        df_row$ci_lower <- ci_df$ci_lower
        df_row$ci_upper <- ci_df$ci_upper
      } else {
        df_row$ci_lower <- NA
        df_row$ci_upper <- NA
      }

      rows <- append(rows, list(df_row))
    }
  }

  # Compute outcome shares for treated units at initial period
  outcome_shares <- compute_outcome_shares(model, pmf_type = "initial")

  # Combine and format
  plot_data_df <- do.call(rbind, rows) %>%
    dplyr::mutate(
      j = factor_j_with_priors(j, model$priors),
      initial_outcome = factor_outcome_with_share(initial_outcome, Y_names, outcome_shares),
      type = paste(initial_outcome, j, sep = "_")
    )

  return(plot_data_df)
}

#' Plot Counterfactual Differences by Initial Outcomes (from Period 0)
#'
#' This function plots the difference (treated - control) in expected outcomes
#' conditional on initial outcomes at the first period (t=0). Uses non-Markovian
#' transition matrices Ps_control_from_0 and Ps_treated_from_0.
#'
#' @param model A fitted transition model object with Ps_control_from_0 and Ps_treated_from_0.
#' @param uniform_ci A logical; if TRUE, compute uniform confidence intervals.
#' @param recenter_t A logical; if TRUE, recenter time relative to the treated period.
#' @param include_post_treatment_periods A logical; if FALSE (default), only include
#'   pre-treatment periods and compute uniform CIs only on those periods.
#'
#' @return A ggplot object.
#' @export
plot_counterfactuals_differences_by_initial_outcomes <- function(
    model,
    uniform_ci = TRUE,
    recenter_t = TRUE,
    include_post_treatment_periods = FALSE) {

  # Extract data
  J <- length(model$priors)
  treated_period <- model$treated_period

  # Generate differences data with CIs
  plot_data_df <- generate_plot_data_df_for_counterfactuals_differences_initial_outcomes(
    model, uniform_ci = uniform_ci,
    include_post_treatment_periods = include_post_treatment_periods
  )

  # Recenter time if requested (relative to treated period)
  # 0 = last pre-treatment period, 1 = first post-treatment period
  if (recenter_t) {
    plot_data_df <- plot_data_df %>%
      mutate(t = t - (treated_period - 1))
  }

  # Check if CIs are available
  has_ci <- !all(is.na(plot_data_df$ci_lower))

  # Create interaction for dodging (combination of initial_outcome and j)
  plot_data_df <- plot_data_df %>%
    mutate(dodge_group = interaction(initial_outcome, j))

  # Get colors from theme
  linecolors <- get_line_colors(unique(plot_data_df$j))

  # Position dodge for spacing between points at same time
  pd <- ggplot2::position_dodge(width = 0.3)

  # Build plot - color by j, shape by initial_outcome
  p <- ggplot2::ggplot(plot_data_df, ggplot2::aes(
    x = t,
    y = estimate,
    color = j,
    shape = initial_outcome,
    group = dodge_group
  )) +
    ggplot2::geom_point(size = 3, position = pd)

  # Add CI error bars if available
  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2,
      position = pd
    )
  }

  # Add reference lines and styling
  p <- p +
    reference_hline() +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(values = linecolors) +
    ggplot2::labs(
      x = "Year Relative to the Policy Change",
      y = "Difference (Treated - Control)",
      color = "Latent Type",
      shape = "Initial Outcome (t=1)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 1),
      color = if (J == 1) "none" else ggplot2::guide_legend(order = 2)
    )

  # Save plot
  save_plot(
    get_plot_counterfactuals_by_pretreatment_outcomes_directory(
      model$empirical_spec, J, suffix = "-differences-from-initial"
    ),
    p
  )

  # Return with title
  p + ggplot2::ggtitle(
    empirical_spec_to_plot_title(
      model$empirical_spec,
      "Counterfactual Differences by Initial Outcomes (from t=1)"
    )
  )
}

#' Plot Counterfactual Differences Comparison (Initial vs Pretreatment Outcomes)
#'
#' Creates a two-panel comparison plot showing counterfactual differences
#' (treated - control) conditioned on:
#' - (i) Initial outcomes (t=1) in the left panel
#' - (ii) Pretreatment outcomes in the right panel
#'
#' Both panels share a common legend at the bottom.
#'
#' @param model A fitted transition model with bootstrap_estimates for CIs.
#' @param uniform_ci Logical; if TRUE, use uniform CIs. Default TRUE.
#' @param recenter_t Logical; if TRUE, recenter time so treatment period is 0.
#' @param include_post_treatment_periods Logical; if TRUE, include post-treatment
#'   periods in the initial outcomes panel. Default FALSE.
#'
#' @return A combined ggplot object with two panels.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom ggpubr ggarrange
#' @importFrom magrittr %>%
#' @export
plot_counterfactuals_comparison_differences <- function(
    model,
    uniform_ci = TRUE,
    recenter_t = TRUE,
    include_post_treatment_periods = FALSE,
    x_lab = NULL) {

  J <- length(model$priors)
  treated_period <- model$treated_period

  # Use model's empirical_spec$x_lab if available and x_lab not provided

  if (is.null(x_lab)) {
    x_lab <- if (!is.null(model$empirical_spec) && !is.null(model$empirical_spec$x_lab)) {
      model$empirical_spec$x_lab
    } else {
      "Year Relative to the Policy Change"
    }
  }

  # --- Panel 1: Counterfactual Differences by Initial Outcomes ---
  # Note: j is already factored with labels inside generate_plot_data_df_for_counterfactuals_differences_initial_outcomes
  plot_data_initial <- generate_plot_data_df_for_counterfactuals_differences_initial_outcomes(
    model, uniform_ci = uniform_ci,
    include_post_treatment_periods = include_post_treatment_periods
  )

  if (recenter_t) {
    plot_data_initial <- plot_data_initial %>%
      mutate(t = t - (treated_period - 1))
  }

  has_ci_initial <- !all(is.na(plot_data_initial$ci_lower))
  plot_data_initial <- plot_data_initial %>%
    mutate(dodge_group = interaction(initial_outcome, j))

  # Get colors and name them according to factor levels
  j_levels <- levels(plot_data_initial$j)
  linecolors <- get_line_colors(j_levels)
  names(linecolors) <- j_levels
  pd <- ggplot2::position_dodge(width = 0.3)

  p_initial <- ggplot2::ggplot(plot_data_initial, ggplot2::aes(
    x = t, y = estimate, color = j, shape = initial_outcome, group = dodge_group
  )) +
    ggplot2::geom_point(size = 3, position = pd)

  if (has_ci_initial) {
    p_initial <- p_initial + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2, position = pd
    )
  } else {
    # Only add lines when no CIs (to avoid clutter)
    p_initial <- p_initial + ggplot2::geom_line(position = pd, alpha = 0.5)
  }

  p_initial <- p_initial +
    reference_hline() +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(values = linecolors) +
    ggplot2::labs(
      title = "(i) Counterfactuals by Initial Outcomes",
      x = x_lab,
      y = "Difference (Treated - Control)",
      color = "Latent Type",
      shape = "Conditioning Outcome"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = ggplot2::unit(-3, "pt"),
      plot.title = ggplot2::element_text(size = 10, hjust = 0.5)
    ) +
    ggplot2::scale_x_continuous(breaks = integer_breaks()) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 1, nrow = 1),
      color = if (J == 1) "none" else ggplot2::guide_legend(order = 2, nrow = 1)
    )

  # --- Panel 2: Counterfactual Differences by Pretreatment Outcomes ---
  # Use Ps_control_empirical for consistency with generate_plot_data_df_for_counterfactuals_differences
  ts_post_treatment <- (treated_period - 1):(length(model$Ps_control_empirical[[1]]))
  ts <- c(treated_period - 1, (ts_post_treatment + 1))

  # Note: j is already factored with labels inside generate_plot_data_df_for_counterfactuals_differences
  plot_data_pre <- generate_plot_data_df_for_counterfactuals_differences(
    model, ts, uniform_ci = uniform_ci
  )

  if (recenter_t) {
    plot_data_pre <- plot_data_pre %>%
      mutate(t = t - (treated_period - 1))
  }

  # Filter out t=0 where difference is 0 by construction
  plot_data_pre <- plot_data_pre %>%
    filter(t >= ifelse(recenter_t, 1, treated_period))

  has_ci_pre <- !all(is.na(plot_data_pre$ci_lower))
  plot_data_pre <- plot_data_pre %>%
    mutate(dodge_group = interaction(initial_outcome, j))

  p_pre <- ggplot2::ggplot(plot_data_pre, ggplot2::aes(
    x = t, y = estimate, color = j, shape = initial_outcome, group = dodge_group
  )) +
    ggplot2::geom_point(size = 3, position = pd)

  if (has_ci_pre) {
    p_pre <- p_pre + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2, position = pd
    )
  } else {
    # Only add lines when no CIs (to avoid clutter)
    p_pre <- p_pre + ggplot2::geom_line(position = pd, alpha = 0.5)
  }

  p_pre <- p_pre +
    reference_hline() +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(values = linecolors) +
    ggplot2::labs(
      title = "(ii) Counterfactuals by Pretreatment Outcomes",
      x = x_lab,
      y = "Difference (Treated - Control)",
      color = "Latent Type",
      shape = "Conditioning Outcome"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = ggplot2::unit(-3, "pt"),
      plot.title = ggplot2::element_text(size = 10, hjust = 0.5)
    ) +
    ggplot2::scale_x_continuous(breaks = integer_breaks()) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 1, nrow = 1),
      color = if (J == 1) "none" else ggplot2::guide_legend(order = 2, nrow = 1)
    )

  # --- Combine panels with common legend ---
  combined_plot <- ggpubr::ggarrange(
    p_initial, p_pre,
    ncol = 2, nrow = 1,
    common.legend = TRUE,
    legend = "bottom"
  )

  # Save plot
  save_plot(
    get_plot_counterfactuals_by_pretreatment_outcomes_directory(
      model$empirical_spec, J, suffix = "-comparison-differences"
    ),
    combined_plot
  )

  combined_plot
}

