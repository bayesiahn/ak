

#' Plot Comparison of Counterfactuals Implied by Models
#'
#' This function generates a plot comparing treated and counterfactual outcomes.
#'
#' @param empirical_analysis_result An object containing the results of the empirical analysis.
#' @param show_did Logical, whether to show DiD counterfactuals.
#' @param show_ours Logical, whether to show our method's counterfactuals.
#' @return A ggplot object showing the comparison.
#' @export
plot_all_counterfactual_comparison <- function(empirical_analysis_result,
                                               show_did = TRUE,
                                               show_ours = TRUE) {
  # Extract transition models and did_model from the empirical analysis result
  models <- empirical_analysis_result$transition_models
  did_model <- empirical_analysis_result$did_result
  did_conditional_lag_models <- empirical_analysis_result$did_conditional_lag_models
  # Exclude lag = 1 model
  did_conditional_lag_models <- did_conditional_lag_models[
    sapply(did_conditional_lag_models, function(did_conditional_lag_model) did_conditional_lag_model$lag != 1)]

  # Extract outcome name and variable name
  empirical_spec <- models[[1]]$empirical_spec
  outcome_name <- empirical_spec$outcome_name
  outcome_varname <- empirical_spec$outcome_varname
  x_lab <- models[[1]]$empirical_spec$x_lab

  y_g <- get_y_numeric_and_g_from_model(models[[1]])
  y_numeric <- y_g$y_numeric
  g <- y_g$g
  treated_period <- models[[1]]$treated_period

  # Get counterfactuals
  display_J <- length(models) > 1
  df_counterfactuals <- lapply(models, get_df_counterfactual_transition_model, y_numeric, g, display_J) %>%
    bind_rows() %>%
    rbind(get_df_counterfactual_did_model(did_model, y_numeric, g)) %>%
    rbind(get_df_counterfactual_did_conditional_lag_models(did_conditional_lag_models, y_numeric, g)) %>%
    rbind(get_df_outcome_treated(y_numeric, g)) %>%
    mutate(status = ifelse(type == "Treated", "Treated", "Control"))  %>%
    filter(t >= treated_period - 1)
  types <- unique(sort(df_counterfactuals$type))
  Js <- sapply(models, function(x) length(x$priors))


  if (!show_did) {
    df_counterfactuals <- df_counterfactuals %>%
      filter(type != TYPE_NAME_COUNTERFACTUAL_DID)
  }
  if (!show_ours) {
    # remove any type that contains TYPE_NAME_COUNTERFACTUAL_OURS
    df_counterfactuals <- df_counterfactuals %>%
      filter(!grepl("Our Method", type)) %>%
      filter(type != paste(TYPE_NAME_COUNTERFACTUAL_OURS, "Our Method, J = 1"))
  }

  # Generate aggregate plots
  p <- df_counterfactuals %>%
    dplyr::mutate(t = as.factor(t - treated_period + 1)) %>%
    ggplot2::ggplot(ggplot2::aes(x = t, y = y, group = type, shape = type, color = type, linetype = type)) +
    ggplot2::geom_line(size = 1.05) +
    ggplot2::geom_point(size = 2.5) +
    ak::get_jtpa_plot_theme(outcome_name) +
    ggplot2::scale_color_manual(values= get_line_colors_counterfactual(types, Js),
                                name = "Model") +
    ggplot2::scale_linetype_manual(values=get_line_types_counterfactual(types, Js),
                                   name = "Model") +
    ggplot2::scale_shape_manual(values = get_point_shapes_counterfactual(types, Js),
                                name = "Model") +
    ggplot2::xlab(x_lab) +
    add_vertical_lines_at_g(1.0) +
    reference_hline()

  p_with_title <- p +
    ggplot2::ggtitle(empirical_spec_to_plot_title(empirical_spec, "Counterfactuals Comparison"))

  save_plot(get_plot_counterfactuals_comparison_directory(empirical_spec), p)
  save_plot(get_plot_counterfactuals_comparison_directory(empirical_spec, "-with-title"), p_with_title)

  p_with_title

}

#' Plot Simplified Comparison of Counterfactuals
#'
#' This function generates a simplified plot comparing treated outcomes and
#' counterfactuals from DiD and Our Model (J=1 only).
#'
#' @param empirical_analysis_result An object containing the results of the empirical analysis.
#' @param suffix Optional suffix for the saved filename.
#' @return A ggplot object showing the comparison.
#' @export
plot_counterfactual_comparison_simple <- function(empirical_analysis_result,
                                                  suffix = "-simple") {
  # Extract transition models and did_model from the empirical analysis result
  models <- empirical_analysis_result$transition_models
  did_model <- empirical_analysis_result$did_result

  # Only keep J=1 model
  J1_model <- models[[which(sapply(models, function(x) length(x$priors)) == 1)[1]]]

  # Extract outcome name and variable name
  empirical_spec <- J1_model$empirical_spec
  outcome_name <- empirical_spec$outcome_name
  outcome_varname <- empirical_spec$outcome_varname
  x_lab <- J1_model$empirical_spec$x_lab

  y_g <- get_y_numeric_and_g_from_model(J1_model)
  y_numeric <- y_g$y_numeric
  g <- y_g$g
  treated_period <- J1_model$treated_period

  # Get counterfactuals with simplified labels
  df_counterfactual_ours <- get_df_counterfactual_transition_model(
    J1_model, y_numeric, g, display_J = FALSE, type_name_counterfactual_ours = "Our Model")
  df_counterfactual_did <- get_df_counterfactual_did_model(did_model, y_numeric, g, "DiD")
  df_treated <- get_df_outcome_treated(y_numeric, g, "Treated")

  df_counterfactuals <- rbind(df_counterfactual_ours, df_counterfactual_did, df_treated) %>%
    mutate(status = ifelse(type == "Treated", "Treated", "Control")) %>%
    filter(t >= treated_period - 1)

  types <- c("Treated", "DiD", "Our Model")

  # Define colors, linetypes and shapes for simplified plot
  colors <- c("Treated" = "black", "DiD" = "red", "Our Model" = "blue")
  linetypes <- c("Treated" = "solid", "DiD" = "dashed", "Our Model" = "dashed")
  shapes <- c("Treated" = 16, "DiD" = 1, "Our Model" = 17)

  # Generate aggregate plots
  p <- df_counterfactuals %>%
    mutate(type = factor(type, levels = types)) %>%
    dplyr::mutate(t = as.factor(t - treated_period + 1)) %>%
    ggplot2::ggplot(ggplot2::aes(x = t, y = y, group = type, shape = type, color = type, linetype = type)) +
    ggplot2::geom_line(size = 1.05) +
    ggplot2::geom_point(size = 2.5) +
    ak::get_jtpa_plot_theme(outcome_name) +
    ggplot2::scale_color_manual(values = colors, name = "Model") +
    ggplot2::scale_linetype_manual(values = linetypes, name = "Model") +
    ggplot2::scale_shape_manual(values = shapes, name = "Model") +
    ggplot2::xlab(x_lab) +
    add_vertical_lines_at_g(1.0) +
    reference_hline()

  p_with_title <- p +
    ggplot2::ggtitle(empirical_spec_to_plot_title(empirical_spec, "Counterfactuals Comparison"))

  save_plot(get_plot_counterfactuals_comparison_directory(empirical_spec, suffix), p)
  save_plot(get_plot_counterfactuals_comparison_directory(empirical_spec, paste0(suffix, "-with-title")), p_with_title)

  p_with_title
}


#' Plot Comparison of Average Treatment Effects (ATT)
#'
#' This function generates a plot comparing ATT implied by different models.
#'
#' @param empirical_analysis_result An object containing the results of the empirical analysis.
#' @param param_name A character string for the parameter name to plot (default is "ATT").
#' @param display_ci Logical, whether to display confidence intervals (default is TRUE).
#' @return A ggplot object showing the comparison of ATT estimates.
#' @export
plot_all_att_comparison <- function(empirical_analysis_result,
                                    param_name = "ATT",
                                    display_ci = TRUE,
                                    suffix = "") {
  # Extract transition models and did_model from the empirical analysis result
  models <- empirical_analysis_result$transition_models
  did_model <- empirical_analysis_result$did_result
  did_conditional_lag_models <- empirical_analysis_result$did_conditional_lag_models
  # Exclude lag = 1 model
  did_conditional_lag_models <- did_conditional_lag_models[
    sapply(did_conditional_lag_models, function(did_conditional_lag_model) did_conditional_lag_model$lag != 1)]

  # Extract outcome name and variable name
  empirical_spec <- empirical_analysis_result$empirical_spec
  outcome_name <- empirical_spec$outcome_name
  outcome_varname <- empirical_spec$outcome_varname
  x_lab <- empirical_spec$x_lab

  # Construct df_att
  Js <- sapply(models, function(x) length(x$priors))
  df_att <- lapply(models, get_df_att_transition_model, param_name) %>%
    do.call(rbind, .) %>%
    rbind(get_df_att_did_model(did_model)) %>%
    rbind(get_df_att_did_conditional_lag_models(did_conditional_lag_models)) %>%
    dplyr::rename(ATT = estimate) %>%
    mutate(t = as.numeric(t))
  has_ci <- !(any(is.na(df_att$ci_lower)) || any(is.na(df_att$ci_upper)))

  types <- unique(df_att$type)

  # Generate aggregate plots with confidence intervals
  suppressWarnings({p <- df_att %>%
    ggplot2::ggplot(ggplot2::aes(x = t, y = !!sym(param_name), shape = type,
                                 group = type, color = type)) +
    ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(0.5)) +
    ak::get_jtpa_plot_theme(outcome_name) +
    ggplot2::labs(y = "Average treatment effects (ATT, aggregated)") +
    ggplot2::scale_color_manual(
      values = get_line_colors_comparison(types, Js),
      name = "Model"
    ) +
    ggplot2::scale_shape_manual(values = get_point_shapes_comparison(types, Js),
                                name = "Model") +
    ggplot2::xlab(x_lab) +
    ggplot2::labs(shape = "Model") +
    ggplot2::scale_x_continuous(breaks = integer_breaks())
  })

  if (has_ci && display_ci) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        position = ggplot2::position_dodge(0.5),
        width = 0.2
      )
  }

  p <- p +
    reference_hline()

  p_with_title <- p +
    ggplot2::ggtitle(empirical_spec_to_plot_title(empirical_spec, "ATT Comparison"))

  save_plot(get_plot_att_comparison_directory(empirical_spec), p)
  save_plot(get_plot_att_comparison_directory(empirical_spec, "-with-title"), p_with_title)

  p_with_title
}

#' Plot Simplified Comparison of Average Treatment Effects (ATT)
#'
#' This function generates a simplified plot comparing ATT from DiD and Our Model (J=1 only).
#'
#' @param empirical_analysis_result An object containing the results of the empirical analysis.
#' @param param_name A character string for the parameter name to plot (default is "ATT").
#' @param display_ci Logical, whether to display confidence intervals (default is TRUE).
#' @param suffix Optional suffix for the saved filename.
#' @return A ggplot object showing the comparison of ATT estimates.
#' @export
plot_att_comparison_simple <- function(empirical_analysis_result,
                                       param_name = "ATT",
                                       display_ci = TRUE,
                                       suffix = "-simple") {
  # Extract transition models (only J=1) and did_model from the empirical analysis result
  models <- empirical_analysis_result$transition_models
  did_model <- empirical_analysis_result$did_result

  # Only keep J=1 model
  J1_model <- models[[which(sapply(models, function(x) length(x$priors)) == 1)[1]]]

  # Extract outcome name and variable name
  empirical_spec <- empirical_analysis_result$empirical_spec
  outcome_name <- empirical_spec$outcome_name
  outcome_varname <- empirical_spec$outcome_varname
  x_lab <- empirical_spec$x_lab

  # Construct df_att with renamed labels
  df_att_j1 <- get_df_att_transition_model(J1_model, param_name) %>%
    mutate(type = "Our Model")
  df_att_did <- get_df_att_did_model(did_model) %>%
    mutate(type = "DiD")

  df_att <- rbind(df_att_j1, df_att_did) %>%
    dplyr::rename(ATT = estimate) %>%
    mutate(t = as.numeric(t))
  has_ci <- !(any(is.na(df_att$ci_lower)) || any(is.na(df_att$ci_upper)))

  types <- c("DiD", "Our Model")

  # Define colors and shapes for simplified plot
  colors <- c("DiD" = "red", "Our Model" = "blue")
  shapes <- c("DiD" = 16, "Our Model" = 17)

  # Generate aggregate plots with confidence intervals
  suppressWarnings({p <- df_att %>%
    mutate(type = factor(type, levels = types)) %>%
    ggplot2::ggplot(ggplot2::aes(x = t, y = !!sym(param_name), shape = type,
                                 group = type, color = type)) +
    ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(0.5)) +
    ak::get_jtpa_plot_theme(outcome_name) +
    ggplot2::labs(y = "Average treatment effects (ATT, aggregated)") +
    ggplot2::scale_color_manual(
      values = colors,
      name = "Model"
    ) +
    ggplot2::scale_shape_manual(values = shapes,
                                name = "Model") +
    ggplot2::xlab(x_lab) +
    ggplot2::labs(shape = "Model") +
    ggplot2::scale_x_continuous(breaks = integer_breaks())
  })

  if (has_ci && display_ci) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        position = ggplot2::position_dodge(0.5),
        width = 0.2
      )
  }

  p <- p +
    reference_hline()

  p_with_title <- p +
    ggplot2::ggtitle(empirical_spec_to_plot_title(empirical_spec, "ATT Comparison"))

  save_plot(get_plot_att_comparison_directory(empirical_spec, suffix), p)
  save_plot(get_plot_att_comparison_directory(empirical_spec, paste0(suffix, "-with-title")), p_with_title)

  p_with_title
}

