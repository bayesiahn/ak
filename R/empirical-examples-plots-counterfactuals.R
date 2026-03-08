# Line color specifications for plots
#' Line Colors for Counterfactual Models
#'
#' These constants define the line colors used for displaying outcomes in plots.
#'
LINECOLORS_DID <- c("red")
LINECOLORS_OURS <- colorRampPalette(c("#1E88E5", "#E69F00"))(6)
LINECOLORS_MORE_LAG <- "#60ccaf"
LINECOLORS_MORE_LAGS <- colorRampPalette(c("#60ccaf", "#000000"))(5)
SHAPETYPES_DID <- 19
SHAPETYPES_OURS <- c(15, 17, 18, 25)
SHAPETYPES_OURS_COUNTERFACTUALS <- c(0, 2, 5, 6, 7, 9, 10, 11, 12, 13, 14, 3, 4)
SHAPETYPES_MORE_LAG_COUNTERFACTUALS <- 14

#' Get Line Colors for Model Comparison Plots
#'
#' Returns a named vector of line colors for comparing different model types
#' in counterfactual plots.
#'
#' @param types Character vector of model type labels.
#' @param Js Integer vector of latent group counts for each model.
#'
#' @return Named character vector of hex color codes.
#' @keywords internal
get_line_colors_comparison <- function(types, Js) {
  linecolors <- c(LINECOLORS_OURS[Js], LINECOLORS_DID)

  if (has_did_conditional_lag_model(types)) {
    # Add color for more lag model
    linecolors <- c(linecolors,
                    LINECOLORS_MORE_LAGS[1:(count_did_conditional_lag_models(types))])
  }

  setNames(linecolors, types)
}

#' Get Point Shapes for Model Comparison Plots
#'
#' Returns a named vector of point shapes for comparing different model types
#' in counterfactual plots.
#'
#' @param types Character vector of model type labels.
#' @param Js Integer vector of latent group counts for each model.
#'
#' @return Named numeric vector of ggplot2 shape codes.
#' @keywords internal
get_point_shapes_comparison <- function(types, Js) {
  linecolors <- c(SHAPETYPES_OURS[Js], SHAPETYPES_DID)

  if (has_did_conditional_lag_model(types)) {
    # Add color for more lag model
    linecolors <- c(linecolors,
                    rep(SHAPETYPES_MORE_LAG_COUNTERFACTUALS, count_did_conditional_lag_models(types)))
  }
  setNames(linecolors, types)
}

#' Generate Line Colors for Counterfactual Types
#'
#' This function assigns colors to each type of counterfactual model for plotting purposes.
#'
#' @param types A character vector of model types.
#' @param Js A numeric vector indicating the number of models.
#' @return A named vector of colors for each model type.
#' @keywords internal
get_line_colors_counterfactual <- function(types, Js) {
  # Generate colors for our method
  linecolors_ours <- c(LINECOLORS_OURS[Js])
  linecolors_ours_after_j1 <- linecolors_ours[which(Js > 1)]
  if (has_did_conditional_lag_model(types)) {
    # Add color for more lag model
    linecolors_ours <- c(linecolors_ours[1], LINECOLORS_MORE_LAGS[1:(count_did_conditional_lag_models(types))], linecolors_ours_after_j1)
  }

  # Combine all colors
  linecolors <- c(LINECOLORS_DID, linecolors_ours, "black")
  # Return named color vector
  setNames(linecolors, types)
}
#' Check if Conditional Lag Models are Present
#'
#' Checks whether any model type labels contain "Lags", indicating the presence
#' of conditional lag DiD models.
#'
#' @param types Character vector of model type labels.
#'
#' @return Logical indicating whether any conditional lag models are present.
#' @keywords internal
has_did_conditional_lag_model <- function(types) {
  count_did_conditional_lag_models(types) > 0
}

#' Count Conditional Lag Models
#'
#' Counts the number of model type labels containing "Lags", which indicate
#' conditional lag DiD models.
#'
#' @param types Character vector of model type labels.
#'
#' @return Integer count of conditional lag models.
#' @keywords internal
count_did_conditional_lag_models <- function(types) {
  sum(grepl("Lags", types, ignore.case = TRUE))
}


#' Generate Line Types for Counterfactual Types
#'
#' This function assigns line types to each type of counterfactual model for plotting.
#'
#' @param types A character vector of model types.
#' @param Js A numeric vector indicating the number of models.
#' @return A named vector of line types for each model type.
#' @keywords internal
get_line_types_counterfactual <- function(types, Js) {
  # Define line types
  linetypes <- rep("dashed", length(Js) + 1)
  if (has_did_conditional_lag_model(types)) {
    # Add line type for more lag model
    linetypes <- c(rep("dashed", count_did_conditional_lag_models(types)), linetypes)
  }
  linetypes <- c(linetypes, "solid")
  # Return named line type vector
  setNames(linetypes, types)
}

#' Generate Point Shapes for Counterfactual Types
#'
#' Assigns point shapes to each type of counterfactual model for plotting.
#'
#' @param types A character vector of model types.
#' @param Js A numeric vector indicating the number of models.
#' @return A named vector of ggplot2 shape codes for each model type.
#' @keywords internal
get_point_shapes_counterfactual <- function(types, Js) {
  # Define line types
  shape_did <- 1
  shape_treated <- 16
  shapes_ours <- rep(SHAPETYPES_OURS_COUNTERFACTUALS, length(Js))[1:length(Js)]
  shapes_ours_after_j1 <- shapes_ours[which(Js > 1)]
  if (has_did_conditional_lag_model(types)) {
    # Add shape for more lag model
    shapes_ours <- c(shapes_ours[1], rep(SHAPETYPES_MORE_LAG_COUNTERFACTUALS, count_did_conditional_lag_models(types)), shapes_ours_after_j1)
  }

  # Return named line type vector
  shapes <- c(shape_did, shapes_ours, shape_treated)
  setNames(shapes, types)
}


#' Prepare DataFrame for Counterfactual Analysis
#'
#' This function creates a DataFrame for plotting counterfactual outcomes from the provided models.
#'
#' @param model The main counterfactual model.
#' @param did_model The DiD counterfactual model.
#' @param y Outcome variable matrix.
#' @param g Group indicator vector.
#' @param display_J Logical, whether to display the number of models.
#' @return A DataFrame of counterfactual outcomes.
#' @keywords internal
get_df_counterfactuals <- function(model, did_model, y, g, display_J = FALSE) {
  # Number of time periods and models
  T_max <- length(model$Ps_control_empirical[[1]]) + 1
  J <- length(model$priors)
  treated_period <- model$treated_period

  # Compute weighted means by group, time, and treatment
  df_outcome <- ak::get_weighted_means_by_gjt(y, g, rep(1, length(g))) %>%
    rename(y = weighted_mean)
  # Separate treated and control outcomes
  df_outcome_treated <- df_outcome %>% filter(g != 0) %>% ungroup() %>% select(t, y)
  df_outcome_control <- df_outcome %>% filter(g == 0) %>% ungroup() %>% select(t, y)

  # Extract ATT values for both models
  model_summary <- summarize_transition_model(model)
  att_dbt <- c(0, sapply(list(model_summary), "[[", "ATT"))
  att_did <- c(0, did_model$att[((length(did_model$att) - length(att_dbt)) + 2):
                                (length(did_model$att))])

  # Time periods for counterfactuals
  ts <- (T_max - length(att_dbt) + 1):T_max

  # Construct counterfactual dataset
  type_name_counterfactual_ours <- paste0("Our Method",
                                          ifelse(display_J, paste0(", J = ", J), ""))
  data.frame(att = c(att_dbt, att_did),
             t = rep(ts, 2),
             type = c(rep(type_name_counterfactual_ours, each = length(ts)),
                      rep(TYPE_NAME_COUNTERFACTUAL_DID, length(ts)))) %>%
    left_join(df_outcome_treated, by = c("t" = "t")) %>%
    mutate(untreated_counterfactual = y - att) %>%
    mutate(y = untreated_counterfactual) %>%
    select(t, y, type) %>%
    rbind(df_outcome_treated %>% mutate(type = TYPE_NAME_TREATED)) %>%
    filter(t >= treated_period - 1)
  shapes <- c(shape_did, shapes_ours, shape_treated)
  setNames(shapes, types)
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
#' plot_empirical_model_counterfactuals_by_pretreatment_outcomes(model)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
plot_empirical_model_counterfactuals_by_pretreatment_outcomes <- function(model) {
  empirical_spec <- model$empirical_spec
  outcome_varname <- empirical_spec$outcome_varname
  outcome_name <- empirical_spec$outcome_name
  x_lab <- empirical_spec$x_lab
  y_lab <- empirical_spec$y_lab_average_1

  p_pretreatment <- plot_counterfactuals_by_pretreatment_outcomes(model) +
    ggplot2::xlab(x_lab) + ggplot2::ylab(y_lab)
  save_plot(get_plot_counterfactuals_by_pretreatment_outcomes_directory(empirical_spec,
                                                                        J = length(model$priors)),
            p_pretreatment)

  p_pretreatment +
    ggplot2::ggtitle(empirical_spec_to_plot_title(
      model$empirical_spec, "Counterfactuals by Pretreatment Outcomes"
    ))
}

