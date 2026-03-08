#' @title Plot Utility Functions
#'
#' @description
#' Helper functions for creating consistent ggplot2 visualizations in ak.
#'
#' @name plot-utils
NULL

#' Color palette for ak plots
#'
#' A colorblind-friendly palette for plot elements.
#' @keywords internal
palette <- c("#0172b2", "#E69F00", "#CC79A7", "#009e73",
             "#56b4e9", "#332288", "#F0E442")

#' Get Line Types for Treatment Groups
#'
#' Returns ggplot2 linetypes based on the number of treatment groups.
#'
#' @param gs A vector of treatment group indicators.
#' @return A character vector of linetypes.
#' @keywords internal
get_line_types <- function(gs) {
  G <- length(unique(gs))
  linetypes <- rep("solid", G)
  if (G == 2) linetypes <- c("dotted", "solid")
  if (G == 3) linetypes <- c("dotted", "dashed", "solid")
  if (G == 4) linetypes <- c("dotted", "dashed", "twodash", "solid")
  if (G > 4 && 0 %in% gs) linetypes <- c("dotted", rep("solid", G-1))
  return (linetypes)
}

#' Get Shape Types for Treatment Groups
#'
#' Returns ggplot2 point shapes based on the number of treatment groups.
#'
#' @param gs A vector of treatment group indicators.
#' @return A numeric vector of shape codes.
#' @keywords internal
get_shape_types <- function(gs) {
  G <- length(unique(gs))
  shapetypes <- rep(19, G)
  if (G == 2) shapetypes <- c(1, 19)
  if (G == 3) shapetypes <- c(1, 19, 17)
  if (G == 4) shapetypes <- c(1, 19, 17, 15)
  if (G > 4 && 0 %in% gs) shapetypes <- c(1, rep(19, G-1))
  return (shapetypes)
}

#' Get Shape Types for Unique Data Points
#'
#' Returns ggplot2 point shapes for distinguishing many unique groups.
#'
#' @param gs A vector of group indicators.
#' @return A numeric vector of shape codes.
#' @keywords internal
get_shape_types_for_unique_data <- function(gs) {
  G <- length(unique(gs))
  shapetypes <- c(15, 17, 18, 19, 16, 21:25, 0:14)
  if (G > length(shapetypes)) {
    shapetypes <- rep(shapetypes, G)
  }
  return (shapetypes[1:G])
}

#' Get Line Colors for Latent Groups
#'
#' Returns colors from the ak palette for latent groups.
#'
#' @param js A vector of latent group indicators.
#' @return A character vector of hex color codes.
#' @keywords internal
get_line_colors <- function(js) {
  J <- length(unique(js))
  linecolors <- palette[1:J]
  if (J > length(palette)) linecolors <- rep(palette, J)[1:J]
  return (linecolors)
}

#' Get Line Sizes for Treatment Groups
#'
#' Returns line widths for plotting.
#'
#' @param gs A vector of treatment group indicators.
#' @return A numeric vector of line sizes.
#' @keywords internal
get_line_sizes <- function(gs) {
  G <- length(unique(gs))
  linesizes <- rep(2, G)
  return (linesizes)
}

#' Get Latent Type Labels with Proportions
#'
#' Creates labels for latent types that include the proportion (prior probability).
#'
#' @param priors A numeric vector of prior probabilities for each latent type.
#' @return A character vector of labels like "1 (85.5%)", "2 (14.5%)", etc.
#' @export
get_latent_type_labels <- function(priors) {
  J <- length(priors)
  labels <- sapply(1:J, function(j) {
    sprintf("%d (%.1f%%)", j, priors[j] * 100)
  })
  return(labels)
}

#' Get Outcome Labels with Shares
#'
#' Creates labels for outcomes that include the share (probability) of each outcome.
#'
#' @param outcome_names A character vector of outcome names.
#' @param outcome_shares A numeric vector of shares/probabilities for each outcome.
#' @return A character vector of labels like "employed (70.5%)", "unemployed (10.2%)", etc.
#' @export
get_outcome_labels <- function(outcome_names, outcome_shares) {
  labels <- sapply(seq_along(outcome_names), function(i) {
    sprintf("%s (%.1f%%)", outcome_names[i], outcome_shares[i] * 100)
  })
  return(labels)
}

#' Create Factor with Outcome Share Labels
#'
#' Converts a character/factor vector of outcome names to a factor with share labels.
#'
#' @param outcomes A character or factor vector of outcome names.
#' @param outcome_names A character vector of all possible outcome names (in order).
#' @param outcome_shares A numeric vector of shares/probabilities for each outcome.
#' @return A factor with levels labeled like "employed (70.5%)", "unemployed (10.2%)", etc.
#' @export
factor_outcome_with_share <- function(outcomes, outcome_names, outcome_shares) {
  # Convert factor to character if necessary
  if (is.factor(outcomes)) {
    # Extract name part from labels like "employed (70.5%)" or just "employed"
    outcomes_char <- as.character(outcomes)
    outcomes <- gsub("\\s*\\(.*", "", outcomes_char)
  }
  labels <- get_outcome_labels(outcome_names, outcome_shares)
  factor(outcomes, levels = outcome_names, labels = labels)
}

#' Extract Outcome Name from Labeled Factor
#'
#' Extracts the outcome name from a labeled factor like "employed (70.5%)".
#'
#' @param outcomes A character or factor vector of outcome names (possibly with labels).
#' @return A character vector with just the outcome names (without percentage labels).
#' @export
extract_outcome_name <- function(outcomes) {
  if (is.factor(outcomes)) {
    outcomes <- as.character(outcomes)
  }
  gsub("\\s*\\(.*", "", outcomes)
}

#' Compute Outcome Shares from Model
#'
#' Computes the marginal share of each outcome for treated units at a given period.
#'
#' @param model A fitted transition model object.
#' @param pmf_type Character string: "pretreatment" or "initial" to specify which PMF to use.
#' @return A numeric vector of outcome shares summing to 1.
#' @export
compute_outcome_shares <- function(model, pmf_type = "pretreatment") {
  J <- length(model$priors)
  priors_treated <- if (!is.null(model$priors_treated)) model$priors_treated else model$priors

  if (pmf_type == "pretreatment") {
    pmfs <- model$pmfs_pre_treatment_treated
  } else if (pmf_type == "initial") {
    pmfs <- model$pmfs_initial_treated
  } else {
    stop("pmf_type must be 'pretreatment' or 'initial'")
  }

  # Compute weighted average of PMFs across latent groups
  outcome_shares <- Reduce(`+`, lapply(seq_len(J), function(j) {
    priors_treated[j] * pmfs[[j]]
  }))

  return(outcome_shares)
}

#' Create Factor with Latent Type Labels
#'
#' Converts a numeric vector of j values to a factor with proportion labels.
#'
#' @param j A numeric or factor vector of latent type indices (1, 2, 3, ...).
#' @param priors A numeric vector of prior probabilities for each latent type.
#' @return A factor with levels labeled like "1 (85.5%)", "2 (14.5%)", etc.
#' @export
factor_j_with_priors <- function(j, priors) {
  # Convert factor to numeric if necessary
  if (is.factor(j)) {
    # Extract numeric part from labels like "1 (100.0%)" or just "1"
    j_char <- as.character(j)
    j <- as.numeric(gsub("\\s.*", "", j_char))
  }
  labels <- get_latent_type_labels(priors)
  factor(j, levels = 1:length(priors), labels = labels)
}

#' Convert Factor to Numeric
#'
#' Safely converts a factor to numeric values. Handles labeled factors
#' like "1 (100.0%)" by extracting the numeric prefix.
#'
#' @param numeric_vector A factor or character vector to convert.
#' @return A numeric vector.
#' @keywords internal
unfactor_numerics <- function(numeric_vector) {
  char_vector <- as.character(numeric_vector)
  # Extract numeric part from labels like "1 (100.0%)" or just "1"
  return(as.numeric(gsub("\\s.*", "", char_vector)))
}

#' Get Factor Levels as Character
#'
#' Extracts factor levels as character strings.
#'
#' @param vector A factor vector.
#' @return A character vector.
#' @keywords internal
unfactor <- function(vector) {
  levels(vector)[vector]
}

#' Add vertical lines for events (g)
#' @export
#' @param g a vector of events
#' @param display whether to display the vertical lines
#' @return a ggplot object
#' @examples
#' g <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' plot_event_study_group(complete_sample, g) + add_vertical_lines_at_g(g)
#' plot_event_study_group(complete_sample, g) + add_vertical_lines_at_g(g, FALSE)
add_vertical_lines_at_g <- function(g, display = TRUE) {
  # if not needed, return NULL
  if (!display)
    return(NULL)

  g_unique <- sort(unique(unfactor_numerics(g)))

  # choose events except first
  g_unique <- g_unique
  if (length(g_unique) > 1) {
    g_unique <- g_unique[-1]
  }
  event_dates <- data.frame(date = g_unique + 0.5)

  # return vertical lines associated with the events
  geom_vline(data = event_dates, aes(xintercept = date), linetype="dashed")
}

#' Add a horizontal reference line
#'
#' Adds a dashed horizontal reference line to a ggplot. Centralizes the style
#' so that all reference lines in the package look consistent.
#'
#' @export
#' @param yintercept the y value for the horizontal line (default 0)
#' @param color line color (default "black")
#' @param linetype line type (default "dashed")
#' @param alpha line transparency (default 0.5)
#' @return a ggplot2 geom_hline layer
#' @examples
#' library(ggplot2)
#' ggplot(data.frame(x = 1:5, y = rnorm(5)), aes(x, y)) +
#'   geom_point() +
#'   reference_hline()
reference_hline <- function(yintercept = 0, color = "black",
                            linetype = "dashed", alpha = 0.5) {
  ggplot2::geom_hline(yintercept = yintercept, color = color,
                      linetype = linetype, alpha = alpha)
}

#' Get a ggplot theme for ak
#' @export
#' @param gs a vector of unique cohorts
#' @param js a vector of unique groups
#' @param shapes a vector of unique shapes
#' @param y_lab the label for the y-axis
#' @return a list of ggplot2::theme objects
#' @examples
#' get_transition_plot_theme(c(0,1,2),c(1,2,3))
get_transition_plot_theme <- function(gs, js, shapes = NULL, y_lab = NULL) {
  gs <- sort(unique(unfactor_numerics(gs)))
  js <- sort(unique((js)))
  y_lab <- ifelse(is.null(y_lab), "Employment probability", y_lab)
  linetype_title <- "Treatment Status"
  cohort_labels <- c("Control", "Treated")
  show_j_labels <- ggplot2::guide_legend(order = 1)
  show_g_labels <- ggplot2::guide_legend(order = 2)
  if (length(gs) == 1 && gs[1] == 0) cohort_labels <- "Control"
  if (length(gs) == 1 && gs[1] > 0) cohort_labels <- "Treated"
  if (length(gs) == 1) {
    linetype_title <- ""
    cohort_labels <- "Difference in Transition Probability"
    show_g_labels <- "none"
  }
  if (length(js) == 1) {
    show_j_labels <- "none"
  }
  linetypes <- get_line_types(gs)
  linecolors <- get_line_colors(js)
  linesizes <- get_line_sizes(gs)

  # Shape manual
  shape_manual <- NULL
  if (is.null(shapes)) {
    # if shapes are not provided, use default shapes
    shapetypes <- get_shape_types(gs)
    shape_labels <- c("Control", "Treated")
    shape_manual <- ggplot2::scale_shape_manual(values=shapetypes, labels=shape_labels)
  }
  # Use nrow = 1 in guides to keep each legend on single row
  show_j_labels_with_nrow <- if (identical(show_j_labels, "none")) "none" else
    ggplot2::guide_legend(order = 1, nrow = 1)
  show_g_labels_with_nrow <- if (identical(show_g_labels, "none")) "none" else
    ggplot2::guide_legend(order = 2, nrow = 1)

  list(ggplot2::theme_bw(),
       ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                      legend.spacing.y = ggplot2::unit(-3, "pt")),
       ggplot2::scale_linetype_manual(values=linetypes, labels=cohort_labels),
       ggplot2::scale_color_manual(values=linecolors),
       ggplot2::scale_size_manual(values=linesizes),
       shape_manual,
       ggplot2::guides(
           shape = show_g_labels_with_nrow,
           linetype = show_g_labels_with_nrow, # Set linetype legend order
           color = show_j_labels_with_nrow,
         ),
       ggplot2::labs(x = "Time", y = y_lab,
                     color = "Latent Type",
                     linetype = linetype_title,
                     shape = linetype_title))

}

#' Get integer breaks for x or y scale
#' @param n the number of breaks
#' @return a function that returns integer breaks
#' @export
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    unique(breaks)
  }
  return(fxn)
}

#' Get a ggplot theme for JTPA set
#' @export
#' @param y_lab the label for the y-axis
#' @return a list of ggplot2::theme objects
#' @examples
#' get_jtpa_plot_theme(c(0,1,2),c(1,2,3))
get_jtpa_plot_theme <- function(y_lab = NULL) {
  y_lab <- ifelse(is.null(y_lab), "Employment probability", y_lab)

  list(ggplot2::theme_bw(),
       ggplot2::theme(legend.position = "bottom"),
       ggplot2::labs(x = "Month", y = y_lab))
}

#' Attach Latent Group Index and Weights to Data Frame
#'
#' Adds latent group index (j) and corresponding weights to each row of a data
#' frame for weighted averaging calculations.
#'
#' @param df Data frame to attach j and weights to.
#' @param j Integer latent group index.
#' @param weights_j Numeric vector of weights for group j.
#' @return Data frame with j and weight columns attached.
#' @keywords internal
#' @examples
#' df <- data.frame(id = c(1,1,1,2,2,2), t = c(1,2,3,1,2,3),
#'                 y = c(1,2,3,4,5,6), g = c(1,1,1,0,0,0))
#' df
#' attach_j_and_weights(df, 1, c(0.6, 0.4))
#' attach_j_and_weights(df, 2, c(0.6, 0.4))
attach_j_and_weights <- function(df, j, weights_j) {
  T <- as.integer(nrow(df) / length(weights_j))
  df %>%
    dplyr::mutate(j = j,
                  weight = rep(weights_j, T))
}


#' Arrange Multiple ggplot2 Plots Horizontally
#'
#' Arranges multiple ggplot2 plot objects side by side in a single row within a single plot area,
#' ensuring they share a common legend at the bottom. This function is useful for comparing related
#' visualizations or displaying complementary information across several plots.
#'
#' @param ... One or more ggplot2 plot objects to be arranged.
#'
#' @return A combined ggplot object with the input plots arranged side by side in a single row
#'         and a common legend at the bottom.
#' @import ggpubr
#' @importFrom ggplot2 guides guide_legend
#' @examples
#' plot1 <- ggplot(mtcars, aes(mpg, disp)) + geom_point(aes(color = factor(cyl))) + theme_minimal()
#' plot2 <- ggplot(mtcars, aes(hp, wt)) + geom_point(aes(color = factor(cyl))) + theme_minimal()
#' plot3 <- ggplot(mtcars, aes(qsec, drat)) + geom_point(aes(color = factor(cyl))) + theme_minimal()
#' combined_plot <- arrange_multiple_plots(plot1, plot2, plot3)
#' print(combined_plot)
#'
#' @export
arrange_plots <- function(...) {
  plots <- list(...)
  # Prepare plots with common guides
  prepared_plots <- lapply(plots, function(plot) {
    plot + guides(nrow = 1, color = guide_legend(order = 1), linetype = guide_legend(order = 2),
                  shape = guide_legend(order = 2))
  })

  # Arrange plots together in a single row
  combined_plot <- ggpubr::ggarrange(plotlist = prepared_plots, ncol = length(plots), nrow = 1,
                                     common.legend = TRUE, legend = "bottom") +
    guides(nrow = 1, color = guide_legend(order = 1), linetype = guide_legend(order = 2),
           shape = guide_legend(order = 2))

  return(combined_plot)
}

#' Attach Confidence Intervals to Plot Data Frame
#'
#' This function attaches confidence intervals (CIs) to a given plot data frame
#' based on a model's bootstrap estimates using the percentile-t (studentized
#' bootstrap) method.
#'
#' The CIs are constructed by:
#' \enumerate{
#'   \item Computing bootstrap standard deviations for each (j, g, t) combination.
#'   \item Computing bootstrap t-statistics: (boot_estimate - estimate) / bootstrap_sd.
#'   \item For uniform CIs: using the max |t| across all periods as the critical value.
#'   \item For pointwise CIs: using the 95th percentile of |t| at each period.
#'   \item Returning symmetric CIs: estimate +/- critical_value * sd.
#' }
#'
#' @param plot_data_df A data frame containing the plot data. Must have columns
#'   `j`, `g`, `t`, and `y` (point estimates).
#' @param model A model object that contains `bootstrap_estimates` used for
#'   calculating confidence intervals.
#' @param parameter_name_control Name of the parameter for the control group.
#' @param parameter_name_treated (Optional) Name of the parameter for the treated
#'   group. Defaults to the value of `parameter_name_control`.
#' @param current_index The index corresponding to the current period (default is 1).
#' @param past_index The index corresponding to the past period (default is 1).
#' @param uniform_ci Logical; if TRUE (default), compute uniform (simultaneous)
#'   confidence intervals; if FALSE, compute pointwise intervals.
#'
#' @return A data frame with confidence intervals (`ci_upper` and `ci_lower`) added.
#'   The CIs are symmetric around the point estimate.
#' @keywords internal
attach_ci_to_plot_data_df <- function(plot_data_df, model,
                                      parameter_name_control,
                                      parameter_name_treated = NULL,
                                      current_index = 1,
                                      past_index = 1,
                                      uniform_ci = TRUE) {
  if (!("bootstrap_estimates" %in% names(model))) {
    return(plot_data_df)
  }

  # Set up parameter names
  if (is.null(parameter_name_treated)) {
    parameter_name_treated <- parameter_name_control
  }

  # Get number of latent groups
  J <- length(model$priors)

  # Parameters that need model summary computation (not directly in model)
  summary_params <- c("LTATTs", "ATT", "Ps_diffs", "Ps_empirical_diffs")
  needs_summary <- parameter_name_control %in% summary_params

  # Helper to extract estimates for a specific (j, g) combination
  extract_estimates_for_jg <- function(m, j, g, param_name, use_summary) {
    if (use_summary) {
      summary <- summarize_transition_model(m)
      fetch_parameter(summary, param_name, j, current_index, past_index)
    } else {
      fetch_parameter(m, param_name, j, current_index, past_index)
    }
  }

  # Build CI data frame for each (j, g) group
  ci_rows <- list()

  for (j in 1:J) {
    # Determine which g values to process for this j
    g_values <- unique(plot_data_df$g[
      plot_data_df$j == j | unfactor_numerics(plot_data_df$j) == j
    ])
    if (length(g_values) == 0) {
      g_values <- unique(plot_data_df$g)
    }

    for (g_val in g_values) {
      # Determine parameter name based on g
      g_numeric <- if (is.factor(g_val)) {
        unfactor_numerics(g_val)
      } else {
        as.numeric(g_val)
      }
      param_name <- if (g_numeric > 0) {
        parameter_name_treated
      } else {
        parameter_name_control
      }

      # Get point estimates for this (j, g)
      point_est <- extract_estimates_for_jg(
        model, j, g_numeric, param_name, needs_summary
      )
      if (is.null(point_est)) next

      # Get bootstrap estimates for this (j, g)
      boot_ests <- lapply(model$bootstrap_estimates, function(bm) {
        extract_estimates_for_jg(bm, j, g_numeric, param_name, needs_summary)
      })

      # Filter out NULL bootstrap estimates
      boot_ests <- Filter(Negate(is.null), boot_ests)
      if (length(boot_ests) == 0) next

      # Use bootstrap_estimates_to_ci_df for t-statistic based CIs
      ci_df <- bootstrap_estimates_to_ci_df(
        point_est, boot_ests, uniform_ci = uniform_ci
      )

      # Get the time periods for this (j, g) from the plot data
      # ci_df$t contains position indices (1, 2, 3, ...)
      # We need to map these to actual time periods in plot_data_df
      j_match <- plot_data_df$j == j | unfactor_numerics(plot_data_df$j) == j
      g_match <- unfactor_numerics(plot_data_df$g) == g_numeric
      relevant_ts <- sort(unique(plot_data_df$t[j_match & g_match]))

      # Map position indices to actual time periods
      n_periods <- length(ci_df$t)
      if (length(relevant_ts) >= n_periods) {
        actual_ts <- relevant_ts[1:n_periods]
      } else {
        # Fallback: use position indices + offset
        actual_ts <- as.numeric(ci_df$t)
      }

      # Create row with j, g, t, ci_lower, ci_upper
      ci_row <- data.frame(
        j = j,
        g = g_numeric,
        t = actual_ts,
        ci_lower = ci_df$ci_lower,
        ci_upper = ci_df$ci_upper
      )

      ci_rows <- append(ci_rows, list(ci_row))
    }
  }

  if (length(ci_rows) == 0) {
    return(plot_data_df)
  }

  # Combine all CI rows
  ci_df_all <- do.call(rbind, ci_rows)

  # Convert j in plot_data_df to numeric for joining
  plot_data_df_join <- plot_data_df %>%
    dplyr::mutate(
      j_numeric = unfactor_numerics(j),
      g_numeric = unfactor_numerics(g),
      t_numeric = as.numeric(t)
    )

  # Join CIs to plot data
  ci_df_all <- ci_df_all %>%
    dplyr::rename(j_numeric = j, g_numeric = g, t_numeric = t)

  plot_data_df <- plot_data_df_join %>%
    dplyr::left_join(ci_df_all, by = c("j_numeric", "g_numeric", "t_numeric")) %>%
    dplyr::select(-j_numeric, -g_numeric, -t_numeric)

  return(plot_data_df)
}

#' Fetch Parameter Values from a Model
#'
#' This function extracts parameter values from a model object. For matrix-based parameters, it retrieves elements using the specified indices.
#'
#' @param model A model object containing the parameter to be fetched.
#' @param parameter_name The name of the parameter to fetch.
#' @param j The group index for which the parameter is fetched (default is 1).
#' @param current_index The index corresponding to the current period (default is 1).
#' @param past_index The index corresponding to the past period (default is 1).
#'
#' @return The fetched parameter value(s) from the model. Returns `NULL` if the parameter does not exist.
#'
#' @examples
#' model <- list(parameter_name = list(list(matrix(c(1, 2, 3, 4), nrow = 2))))
#' fetch_parameter(model, "parameter_name", j = 1, current_index = 1, past_index = 2)
#' @keywords internal
fetch_parameter <- function(model, parameter_name, j = 1,
                            current_index = 1,
                            past_index = 1) {
  if (parameter_name %in% names(model)) {
    fetched <- model[parameter_name][[1]][[j]]
    # If each element is a matrix, then use indices to collect element from each matrix
    if (is.list(fetched) && all(sapply(fetched, is.matrix))) {
      fetched <- sapply(fetched, function(x) x[past_index, current_index])
    }
    return(fetched)
  } else {
    return(NULL)
  }
}

#' Convert Quantiles to Confidence Interval Data Frame
#'
#' This function converts quantile estimates from a model into a data frame containing confidence intervals
#' (CIs) for treated and control groups across different time periods. It uses parameter values extracted
#' from the quantiles object to compute the confidence intervals.
#'
#' @param quantiles A list containing quantile estimates, such as bootstrap quantiles, along with group
#'        information (`priors`, `treated_period`, `Ps_treated`, etc.).
#' @param parameter_name_control The name of the parameter for the control group.
#' @param parameter_name_treated (Optional) The name of the parameter for the treated group. Defaults to
#'        `parameter_name_control` if not specified.
#' @param current_index The index corresponding to the current period in the matrices (default is 1).
#' @param past_index The index corresponding to the past period in the matrices (default is 1).
#'
#' @return A data frame with columns `j` (group index), `g` (treatment indicator), `t` (time period),
#'         and `ci` (confidence interval values).
#'
#' @examples
#' quantiles <- list(
#'   priors = list(1, 2),
#'   treated_period = 1,
#'   Ps_treated = list(matrix(0.1, 5, 5), matrix(0.2, 5, 5)),
#'   control_param = list(list(rep(0.5, 5))),
#'   treated_param = list(list(rep(0.6, 5)))
#' )
#' quantiles_to_ci_df(quantiles, "control_param", "treated_param")
#' @keywords internal
quantiles_to_ci_df <- function(quantiles,
                      parameter_name_control,
                      parameter_name_treated = NULL,
                      current_index = 1,
                      past_index = 1) {
  treated_and_control_together <- FALSE
  if (is.null(parameter_name_treated)) {
    parameter_name_treated <- parameter_name_control
    treated_and_control_together <- TRUE
  }


  J <- length(quantiles$priors)
  treated_period <- quantiles$treated_period
  T_max <- (length(quantiles$Ps_treated[[1]]) + 1)
  param_length <- length(quantiles[parameter_name_control][[1]][[1]])
  first_period <- ifelse(param_length == T_max, 1, 2)

  quantiles_to_ci_df_by_j <- function(j, g) {
    if (g > 0) {
      ci <- fetch_parameter(quantiles, parameter_name_treated, j, current_index, past_index)
    } else {
      ci <- fetch_parameter(quantiles, parameter_name_control, j, current_index, past_index)
    }
    data.frame(j = j, g = g,
               ci = ci)
  }

  df_ci_treated <- do.call(rbind, lapply(1:J, quantiles_to_ci_df_by_j, treated_period))
  df_ci_control <- do.call(rbind, lapply(1:J, quantiles_to_ci_df_by_j, 0))

  if (treated_and_control_together) {
    return(df_ci_treated)
  }
  return(rbind(df_ci_treated, df_ci_control))
}
