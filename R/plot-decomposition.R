#' Convert Inflow Model to Plot Data
#'
#' This function extracts relevant data from an inflow model and formats it for plotting.
#'
#' @param inflow_model A model object representing the inflow process.
#' @param pmfs_per_j A list of per-type PMF vectors, where pmfs_per_j[[j]] gives
#'   Pr(Y_{T0}=y|D=1,Z=j) for each outcome y.
#' @return A data frame formatted for plotting inflow data.
inflow_model_to_plot_data_df <- function(inflow_model, pmfs_per_j) {
  # Extract data
  inflow_outcome <- inflow_model$Y$names[which(inflow_model$Y$values != 0)]
  treated_period <- inflow_model$treated_period
  # Start from treated_period (not treated_period - 1) for actual treatment effects
  ts <- treated_period:(length(inflow_model$Ps_control[[1]]) + 1)

  plot_data <- generate_plot_data_df_for_ltatt_by_pretreatment_outcomes(inflow_model, ts) %>%
    # Extract raw outcome name from labeled factor for matching and filtering
    mutate(initial_outcome_raw = extract_outcome_name(initial_outcome)) %>%
    mutate(outcome = inflow_outcome, color = as.character(initial_outcome_raw),
           linetype = "Inflow", g = treated_period,
           type = paste0(initial_outcome_raw, "_to_", inflow_outcome)) %>%
    # Use type-specific PMF: extract integer j and index into pmfs_per_j
    rowwise() %>%
    mutate(j_int = as.integer(gsub("\\s.*", "", as.character(j))),
           weights = pmfs_per_j[[j_int]][match(initial_outcome_raw, inflow_model$Y$names)]) %>%
    ungroup() %>%
    filter(initial_outcome_raw != inflow_outcome) %>%
    mutate(catt = y, catt_weighted = catt * weights)

  # Add reference points at t = treated_period - 1 with catt = 0
  reference_points <- plot_data %>%
    dplyr::select(initial_outcome, j, outcome, color, linetype, g, type, weights) %>%
    dplyr::distinct() %>%
    dplyr::mutate(t = treated_period - 1, y = 0, catt = 0, catt_weighted = 0)

  dplyr::bind_rows(reference_points, plot_data) %>%
    dplyr::arrange(type, t)
}

#' Convert Outflow Model to Plot Data
#'
#' This function extracts relevant data from an outflow model and formats it for plotting.
#'
#' @param outflow_model A model object representing the outflow process.
#' @param pmfs_per_j A list of per-type PMF vectors, where pmfs_per_j[[j]] gives
#'   Pr(Y_{T0}=y|D=1,Z=j) for each outcome y.
#' @param inflow_model The corresponding inflow model for cross-referencing outcomes.
#' @return A data frame formatted for plotting outflow data.
outflow_model_to_plot_data_df <- function(outflow_model, pmfs_per_j, inflow_model) {
  # Extract data
  inflow_outcome <- inflow_model$Y$names[which(inflow_model$Y$values != 0)]
  outflow_outcome <- outflow_model$Y$names[which(outflow_model$Y$values != 0)]
  treated_period <- outflow_model$treated_period
  # Start from treated_period (not treated_period - 1) for actual treatment effects
  ts <- treated_period:(length(outflow_model$Ps_control[[1]]) + 1)

  plot_data <- generate_plot_data_df_for_ltatt_by_pretreatment_outcomes(outflow_model, ts) %>%
    # Extract raw outcome name from labeled factor for matching and filtering
    mutate(initial_outcome_raw = extract_outcome_name(initial_outcome)) %>%
    mutate(outcome = outflow_outcome, color = as.character(outflow_outcome),
           linetype = "Outflow", g = treated_period,
           type = paste0(initial_outcome_raw, "_to_", outflow_outcome)) %>%
    # Use type-specific PMF: extract integer j and index into pmfs_per_j
    rowwise() %>%
    mutate(j_int = as.integer(gsub("\\s.*", "", as.character(j))),
           weights = pmfs_per_j[[j_int]][match(initial_outcome_raw, inflow_model$Y$names)]) %>%
    ungroup() %>%
    filter(initial_outcome_raw == inflow_outcome) %>%
    mutate(catt = y, catt_weighted = catt * weights)

  # Add reference points at t = treated_period - 1 with catt = 0
  reference_points <- plot_data %>%
    dplyr::select(initial_outcome, j, outcome, color, linetype, g, type, weights) %>%
    dplyr::distinct() %>%
    dplyr::mutate(t = treated_period - 1, y = 0, catt = 0, catt_weighted = 0)

  dplyr::bind_rows(reference_points, plot_data) %>%
    dplyr::arrange(type, t)
}

#' Generate a Decomposition Plot
#'
#' This function generates a decomposition plot using given plot data.
#'
#' @param plot_data_df A data frame containing the plot data.
#' @param y_to_use The variable to use for the y-axis, default is "catt".
#' @param outcome_label_function Optional function to format outcome labels in
#'   the legend. Takes a character vector of outcome names and returns formatted
#'   labels. If NULL, raw outcome names are used.
#' @return A ggplot2 object representing the decomposition plot.
plot_decomposition <- function(plot_data_df, y_to_use = "catt",
                               outcome_label_function = NULL) {
  recenter_t = TRUE
  g_max <- max(plot_data_df$g)
  point_size <- 3

  # Apply outcome label formatting if provided
  if (!is.null(outcome_label_function)) {
    plot_data_df <- plot_data_df %>%
      dplyr::mutate(color = outcome_label_function(color))
  }

  plot_data_df %>%
    dplyr::mutate(g = g, j = as.factor(j)) %>%
    dplyr::mutate(t = as.integer(t - ifelse(recenter_t, g_max - 1, 1))) %>%
    ggplot2::ggplot(ggplot2::aes_string(x = "t", y = y_to_use,
                                        group = "type",
                                        color = "color", linetype = "linetype", linesize = "g",
                                        shape = "color")) +
    ggplot2::geom_point(size = ifelse(is.null(point_size), 2, point_size)) +
    ggplot2::geom_line() +
    reference_hline() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    add_vertical_lines_at_g((plot_data_df$g) - 1 - ifelse(recenter_t, (g_max - 1), 0), display = TRUE) +
    ggplot2::scale_x_continuous(breaks = integer_breaks()) +
    ggplot2::labs(x = "Month Relative to the ACA Introduction",
                  color = "Outcome",
                  linetype = "Type",
                  shape = "Outcome")
}

#' Generate Decomposition Plot Data
#'
#' This function generates data for decomposition plots.
#'
#' @param model A model object.
#' @param outcome_of_interest The primary outcome to focus on, defaults to the last outcome.
#' @return A data frame containing decomposition plot data.
#' @export
model_to_decomposition_plot_data_df <- function(model, outcome_of_interest = NULL) {
  if (is.null(outcome_of_interest)) {
    outcome_of_interest <- model$Y$names[length(model$Y$names)]
  }

  inflow_model <- modify_outcomes_of_interest(model, outcome_of_interest)
  other_outcomes <- setdiff(model$Y$names, outcome_of_interest)
  outflow_models <- lapply(other_outcomes, function(outcome) modify_outcomes_of_interest(model, outcome))

  J <- length(inflow_model$priors)
  # Per-type PMFs: pmfs_per_j[[j]] = Pr(Y_{T0}=y|D=1,Z=j)
  pmfs_per_j <- lapply(1:J, function(j) inflow_model$pmfs_pre_treatment_treated[[j]])

  inflow_model_df <- inflow_model_to_plot_data_df(inflow_model, pmfs_per_j)
  outflow_models_df <- lapply(outflow_models, outflow_model_to_plot_data_df, pmfs_per_j, inflow_model)
  plot_data_df <- dplyr::bind_rows(inflow_model_df, outflow_models_df) %>%
    dplyr::select(-dplyr::any_of(c("j_int", "initial_outcome_raw")))
}

#' Plot Conditional ATT
#'
#' This function generates a plot for conditional average treatment effects (ATT).
#'
#' @param model A model object.
#' @param outcome_of_interest The outcome variable of interest.
#' @param outcome_label_function Optional function to format outcome labels in
#'   the legend. Takes a character vector of outcome names and returns formatted
#'   labels. If NULL, raw outcome names are used.
#' @return A ggplot2 object representing the conditional ATT plot.
#' @export
plot_decomposition_catt <- function(model, outcome_of_interest = NULL,
                                    outcome_label_function = NULL) {
  decomposition_plot_data_df <- model_to_decomposition_plot_data_df(model, outcome_of_interest)
  plot_decomposition(decomposition_plot_data_df, y_to_use = "catt",
                     outcome_label_function = outcome_label_function) +
    ggplot2::ylab(paste("Conditional ATT on", outcome_of_interest))
}

#' Plot Conditional ATT (Weighted)
#'
#' This function generates a plot for weighted conditional average treatment effects (ATT).
#'
#' @param model A model object.
#' @param outcome_of_interest The outcome variable of interest.
#' @param outcome_label_function Optional function to format outcome labels in
#'   the legend. Takes a character vector of outcome names and returns formatted
#'   labels. If NULL, raw outcome names are used.
#' @return A ggplot2 object representing the weighted conditional ATT plot.
#' @export
plot_decomposition_catt_weighted <- function(model, outcome_of_interest = NULL,
                                             outcome_label_function = NULL) {
  decomposition_plot_data_df <- model_to_decomposition_plot_data_df(model, outcome_of_interest)
  plot_decomposition(decomposition_plot_data_df, y_to_use = "catt_weighted",
                     outcome_label_function = outcome_label_function) +
    ggplot2::ylab(paste("Conditional ATT on", outcome_of_interest, "(Weighted)"))
}

#' Plot Decomposition by Latent Type
#'
#' Generates a faceted decomposition plot with one panel per latent type,
#' showing type-specific weighted CATTs decomposed into inflow and outflow
#' components. Each panel uses the type-specific pretreatment PMF
#' Pr(Y_{T0}=y|D=1,Z=j) as weights.
#'
#' @param model A TransitionModel object.
#' @param outcome_of_interest The outcome variable of interest.
#' @param y_to_use The variable to use for the y-axis, default is "catt_weighted".
#' @param outcome_label_function Optional function to format outcome labels in
#'   the legend. Takes a character vector of outcome names and returns formatted
#'   labels. If NULL, raw outcome names are used.
#' @return A ggplot2 object with faceted panels for each latent type.
#' @export
plot_decomposition_by_type <- function(model, outcome_of_interest = NULL,
                                       y_to_use = "catt_weighted",
                                       outcome_label_function = NULL) {
  decomposition_plot_data_df <- model_to_decomposition_plot_data_df(model, outcome_of_interest)
  # Factor j with prior labels for panel titles
  decomposition_plot_data_df <- decomposition_plot_data_df %>%
    dplyr::mutate(j = factor_j_with_priors(j, model$priors))
  plot_decomposition(decomposition_plot_data_df, y_to_use = y_to_use,
                     outcome_label_function = outcome_label_function) +
    ggplot2::facet_wrap(~j) +
    ggplot2::ylab(paste("Conditional ATT on", outcome_of_interest,
                        if (y_to_use == "catt_weighted") "(Weighted)" else ""))
}
