#' Plot group averages
#' @export
#' @param y An N by T matrix of outcomes; each row is unit outcome
#' @param g N-length vector of treatment status; treated if g > 0
#' @param weights An N by J matrix of individual weights of group membership; if NULL assumed to be a single group
#' @param difference_by_treatment Boolean indicating if differences in means by treatment should be computed
#' @param display_ci whether to show confidence intervals (for J = 1)
#' @param display_vertical_event_lines whether to display vertical lines for events
#' @return a ggplot object
#' @export
plot_transition_averages <- function(y, g, weights = NULL,
                               difference_by_treatment = TRUE,
                               display_ci = TRUE,
                               display_vertical_event_lines = TRUE) {
  # Get weighted averages df
  weighted_means_by_gjt <- get_weighted_means_by_gjt(y, g, weights,
                                                     difference_by_treatment = difference_by_treatment)

  # Plot them
  plot_weighted_averages(weighted_means_by_gjt,
                         display_ci = display_ci,
                         display_vertical_event_lines = display_vertical_event_lines,
                         difference_by_treatment = difference_by_treatment)
}


#' Plot group averages
#' @export
#' @param y An N by T matrix of outcomes; each row is unit outcome
#' @param g N-length vector of treatment status; treated if g > 0
#' @param weights An N by J matrix of individual weights of group membership; if NULL assumed to be a single group
#' @param y_name_to A scalar value of y to compute transition probabilities
#' @param y_name_from A scalar value of y in the past to compute transition probabilities
#' @param pretreatment_only Boolean indicating if only pretreatment periods should be plotted
#' @param display_ci whether to show confidence intervals (for J = 1)
#' @param uniform_ci A logical; if \code{TRUE}, compute uniform (simultaneous)
#'   confidence intervals; if \code{FALSE}, compute pointwise intervals.
#' @param display_vertical_event_lines whether to display vertical lines for events
#' @param recenter_t whether to recenter the time axis (g is displayed as 0)
#' @return a ggplot object
#' @export
plot_transition_probabilities <- function(model, y_name_to = NULL, y_name_from = NULL,
                                          pretreatment_only = TRUE,
                                          display_ci = TRUE,
                                          uniform_ci = TRUE,
                                          display_vertical_event_lines = TRUE,
                                          recenter_t = TRUE) {

  # Get estimates df
  Ps_control_df <- summarize_Ps(model$Ps_control_empirical, model$Y$names)
  Ps_treated_df <- summarize_Ps(model$Ps_treated_empirical, model$Y$names)
  Ps_diff_df <- Ps_treated_df %>% mutate(prob = prob - Ps_control_df$prob)

  # Get estimates df for bootstrap models
  bootstrap_models <- model$bootstrap_estimates
  treated_period <- model$treated_period
  has_ci <- FALSE
  # Pre-treatment transitions are from t=1 to t=treated_period-2
  # because t indexes transitions (t=k means transition from period k to k+1),
  # and transition t=treated_period-1 goes INTO the treatment period
  t_until <- ifelse(pretreatment_only, treated_period - 2, Inf)

  if (!is.null(bootstrap_models)) {
    J <- length(model$priors)
    # Filter out bootstrap estimates where empirical Ps are NULL (failed replicates)
    valid_boot <- !sapply(model$bootstrap_estimates, function(b) is.null(b$Ps_control_empirical))
    boot_estimates_valid <- model$bootstrap_estimates[valid_boot]
    Ps_control_bootstrap <- lapply(boot_estimates_valid, "[[", "Ps_control_empirical")
    Ps_treated_bootstrap <- lapply(boot_estimates_valid, "[[", "Ps_treated_empirical")
    Ps_control_bootstrap_df <- lapply(Ps_control_bootstrap, summarize_Ps, model$Y$names)
    Ps_treated_bootstrap_df <- lapply(Ps_treated_bootstrap, summarize_Ps, model$Y$names)
    Ps_diff_bootstrap_df <- mapply(function(treated, control) {
      treated %>% mutate(prob = prob - control$prob)
    }, Ps_treated_bootstrap_df, Ps_control_bootstrap_df, SIMPLIFY = FALSE)

    # For each j, compute uniform ci across all (t, from, to) and save
    df_params_list <- vector("list", J)
    for (j in 1:J) {
      # Filter Ps_diff_df for this specific j
      Ps_diff_df_j <- Ps_diff_df %>% filter(j == !!j)

      # Get indices for the transition of interest (per j)
      estimate_indices_j <- Ps_diff_df_j %>% mutate(index = row_number()) %>%
        filter(from == y_name_from & to == y_name_to & t <= t_until) %>%
        pull(index)

      estimate <- Ps_diff_df_j %>%
        filter(t <= t_until) %>% pull(prob)
      bootstrap_estimates <- lapply(Ps_diff_bootstrap_df, function(df) df %>% filter(j == !!j) %>%
                                      filter(t <= t_until) %>% pull(prob))

      # Compute uniform ci across all (t, from, to) and save
      df_params <- bootstrap_estimates_to_ci_df(estimate, bootstrap_estimates,
                                                uniform_ci = uniform_ci) %>%
        filter(row_number() %in% estimate_indices_j) %>%
        mutate(t = Ps_diff_df_j$t[estimate_indices_j])
      df_params_list[[j]] <- df_params %>% mutate(j = j)
    }
    df_params <- do.call(rbind, df_params_list)
    has_ci <- TRUE
  } else {
    # No bootstrap: create df_params from point estimates only
    J <- length(model$priors)
    df_params_list <- vector("list", J)
    for (j in 1:J) {
      Ps_diff_df_j <- Ps_diff_df %>% filter(j == !!j)
      estimate_indices_j <- Ps_diff_df_j %>% mutate(index = row_number()) %>%
        filter(from == y_name_from & to == y_name_to & t <= t_until) %>%
        pull(index)
      df_params_list[[j]] <- tibble::tibble(
        estimate = Ps_diff_df_j$prob[estimate_indices_j],
        t = Ps_diff_df_j$t[estimate_indices_j],
        j = j
      )
    }
    df_params <- do.call(rbind, df_params_list)
  }
  # For pretreatment_only plots, set g so that the last included transition (t_until)
  # displays at t=0 after recentering. The recentering formula is t - (g - 1),
  # so we need g = t_until + 1 to make t_until map to 0.
  g_for_plot <- ifelse(pretreatment_only, t_until + 1, treated_period)
  plot_data_df <- df_params %>%
    rename(weighted_mean = estimate) %>%
    mutate(t = as.numeric(t), g = as.factor(g_for_plot),
           j = factor_j_with_priors(j, model$priors))



  # if (difference_by_treatment)
  suppressMessages(
  p <- plot_weighted_averages(plot_data_df,
                         display_ci = has_ci,
                         display_vertical_event_lines = display_vertical_event_lines,
                         recenter_t = recenter_t) +
    ggplot2::ylab(paste0(y_name_to,
                 " probability conditional on past (",
                 y_name_from, ")")) +
    scale_shape_discrete(
      labels = function(x) rep("Difference (Treated - Control)", length(x))
    ))
  p
}

#' Plot a data frame of weighted averages
#' @param plot_data_df a data frame with columns of g, j, t, weighted mean, and weighted trend
#' @param difference_by_treatment Boolean indicating if differences in means by treatment should be computed
#' @param display_ci whether to show confidence intervals
#' @param display_vertical_event_lines whether to display vertical lines for events
#' @param recenter_t whether to recenter the time axis (g is displayed as 0)
#' @return a ggplot object
#' @keywords internal
plot_weighted_averages <- function(plot_data_df,
                                   difference_by_treatment = FALSE,
                                   display_ci = FALSE,
                                   display_vertical_event_lines = TRUE,
                                   recenter_t = TRUE) {
  # Generate a plot
  g_max <- max(unfactor_numerics(plot_data_df$g))
  p <- plot_data_df %>%
    dplyr::mutate(t = as.integer(t-ifelse(recenter_t, g_max-1, 1))) %>% # To account for t = 0
    ggplot2::ggplot(ggplot2::aes_string(x = "t", y = "weighted_mean",
                                        group = "interaction(g, j)",
                                        color = "j", linetype = "g", linesize = "g",
                                        shape = "g")) +
    get_transition_plot_theme(plot_data_df$g, plot_data_df$j) +
    add_vertical_lines_at_g(unfactor_numerics(plot_data_df$g) - 1 - ifelse(recenter_t, (g_max-1), 0),
                            display = display_vertical_event_lines) +
    reference_hline() +
    ggplot2::scale_x_continuous(breaks = integer_breaks())

  if (display_ci) {
    # Add confidence intervals; dodge space is used to avoid overlap
    dodge_space <- ifelse(length(unique(plot_data_df$g)) > 0, 0.5, 0)
    if ("ci_lower" %in% names(plot_data_df)) {
      p <- p +
        ggplot2::geom_point(position = ggplot2::position_dodge(dodge_space)) +
        ggplot2::geom_errorbar(linetype = "solid",
                               width= 0.5,
                               position = ggplot2::position_dodge(dodge_space),
                               ggplot2::aes(ymin = ci_lower,
                                            ymax = ci_upper))


    } else {
      p <- p +
        ggplot2::geom_point(position = ggplot2::position_dodge(dodge_space)) +
        ggplot2::geom_errorbar(linetype = "solid",
                              width= 0.5,
                              position = ggplot2::position_dodge(dodge_space),
                              ggplot2::aes(ymin = weighted_mean - 1.96*se,
                                            ymax = weighted_mean + 1.96*se))
    }
  } else {
    p <- p +
      ggplot2::geom_point() +
      ggplot2::geom_line()
  }

  p
}
#' @title compute weighted mean and trend by g, j, and t
#' @param y: N by T wide matrix of outcomes
#' @param g: N-vector of treatment timing for each unit
#' @param weights: N by K matrix of weights; if not provided, set to be 1
#' @param difference_by_treatment: Boolean indicating if differences in means by treatment should be computed
#' @return dataframe with columns of g, j, t, weighted mean, se and weighted trend; se returned NA if J > 1
#' @export
#' @examples
#' y <- matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
#' g <- c(1,0)
#' weights <- matrix(c(0.6, 0.4), nrow = 2, ncol = 2)
#' get_weighted_mean_by_gjt(y, g, weights)
get_weighted_means_by_gjt <- function(y, g, weights = NULL,
                                      difference_by_treatment = FALSE) {
  # If NULL weights are given, replace by 1
  if (is.null(weights)) {
    weights <- matrix(1, nrow = nrow(y), ncol = 1)
  }

  # convert wide df to long df
  long_df <- wide_weighted_sample_to_long_df(y,g,weights) %>%
    dplyr::group_by(id, j)

  # compute weighted means and trends
  N <- nrow(y)
  J <- ifelse(is.null(weights) || !is.matrix(weights), 1, ncol(weights))
  weighted_mean_by_gjt <- long_df %>%
    dplyr::group_by(j, g, t) %>%
    dplyr::summarize(weighted_mean = sum(y*weight)/sum(weight),
                     N_gt = dplyr::n(),
                     var_y = var(y)) %>%
    dplyr::mutate(se = sqrt(var_y/N_gt))

  # If J > 1, do not compute se by default (need bootstrap)
  if (J > 1) {
    weighted_mean_by_gjt$se <- NA
  }

  if (difference_by_treatment) {
    return(get_differences_in_weighted_means_by_g(weighted_mean_by_gjt))
  }

  weighted_mean_by_gjt
}

#' @title compute weighted mean and trend by g, j, and t
#' @param y N by T wide matrix of outcomes
#' @param g N-vector of treatment timing for each unit
#' @param weights N by K matrix of weights
#' @param y_value A scalar value of y to compute transition probabilities
#' @param y_value_past A scalar value of y in the past to compute transition probabilities
#' @param difference_by_treatment: Boolean indicating if differences in means by treatment should be computed
#' @return dataframe with columns of g, j, t, weighted mean, se; NaN if no past value of y_value_past observed and se returned NA if J > 1
#' @examples
#' y <- matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
#' g <- c(1,0)
#' weights <- matrix(c(0.6, 0.4), nrow = 2, ncol = 2)
#' get_weighted_mean_by_gjt(y, g, weights)
get_weighted_transition_probabilities_by_gjt <- function(y, g, weights, y_value, y_value_past,
                                                         difference_by_treatment = FALSE) {
  # convert wide df to long df
  long_df <- wide_weighted_sample_to_long_df(y,g,weights)

  # Construct past values and filter so that only those with past values are included
  long_df <- long_df %>%
    dplyr::group_by(id, j) %>%
    dplyr::mutate(y_past = dplyr::lag(y, n = 1)) %>%
    dplyr::mutate(y = (y %in% y_value)) %>%
    dplyr::mutate(y_conditioned = (y_past %in% y_value_past))

  # compute weighted means and trends
  N <- nrow(y)
  J <- ifelse(is.null(weights) || !is.matrix(weights), 1, ncol(weights))
  weighted_mean_by_gjt <- long_df %>%
    dplyr::group_by(j, g, t) %>%
    dplyr::summarize(weighted_mean = sum(y_conditioned*y*weight)/sum(y_conditioned*weight),
                     N_gt = dplyr::n(),
                     var_y = weighted_mean*(1-weighted_mean),
                     se = sqrt(var_y/N_gt))

  if (difference_by_treatment) {
    return(get_differences_in_weighted_means_by_g(weighted_mean_by_gjt))
  }

  weighted_mean_by_gjt
}

#' @title Compute difference in weighted means by g within each t
#' @param weighted_mean_by_gjt: Dataframe with weighted means by g, j, and t
#' @return Dataframe with columns of j, t, difference in means, and se of the difference
get_differences_in_weighted_means_by_g <- function(weighted_mean_by_gjt) {
  J <- length(unique(weighted_mean_by_gjt$j))

  mean_treated <- weighted_mean_by_gjt %>%
    dplyr::filter(unfactor_numerics(g) > 0) %>%
    dplyr::rename(weighted_mean_treated = weighted_mean,
                  N_treated = N_gt,
                  var_treated = var_y) %>%
    dplyr::select(j, t, weighted_mean_treated, N_treated, var_treated)

  mean_control <- weighted_mean_by_gjt %>%
    dplyr::filter(unfactor_numerics(g) == 0) %>%
    dplyr::rename(weighted_mean_control = weighted_mean,
                  N_control = N_gt,
                  var_control = var_y) %>%
    dplyr::select(j, t, weighted_mean_control, N_control, var_control)

  g_max <- max(unfactor_numerics(mean_treated$g))
  difference_means_by_g <- mean_treated %>%
    dplyr::left_join(mean_control, by = c("j", "t")) %>%
    dplyr::mutate(g = as.factor(g_max),
                  weighted_mean = weighted_mean_treated - weighted_mean_control,
                  se = sqrt((var_treated / N_treated) + (var_control / N_control))) %>%
    dplyr::select(g, j, t, weighted_mean, se)

  return(difference_means_by_g)
}


#' @title convert wide matrix y and g into long dataframe
#' @param y: wide matrix of outcomes
#' @param g: vector of treatment timing for each unit
#' @param weights: matrix of weights
#' @return long dataframe with columns of id, y, g, t, and weight
#' @examples
#' y <- matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
#' g <- c(1,0)
#' weights <- matrix(c(0.6, 0.4), nrow = 2, ncol = 2)
#' wide_weighted_sample_to_long_df(y, g, weights)
wide_weighted_sample_to_long_df <- function(y, g, weights) {
  y <- as.matrix(y)
  g <- as.numeric(g)
  weights <- as.matrix(weights)

  N <- nrow(y)
  T <- ncol(y)
  J <- ncol(weights)

  # generate a long dataframe
  long_df <- cbind(t = 1:T, t(y) %>% as.data.frame() %>% `colnames<-`(1:N)) %>%
    tidyr::pivot_longer(cols = as.character(1:N),
                        names_to = 'id', values_to = 'y') %>%
    dplyr::mutate(g = rep(g, T))

  # attach j index and weights for each j to the df
  do.call(rbind, lapply(1:J, function(j)
    attach_j_and_weights(long_df, j, weights[,j]))) %>%
    dplyr::mutate(j = as.factor(j),
                  g = as.factor(g),
                  t = as.numeric(t))
}

