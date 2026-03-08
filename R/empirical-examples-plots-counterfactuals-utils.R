#' Get Outcome Means for Treated Units by Period
#'
#' Computes the weighted mean of a numeric outcome variable for treated groups (excluding control group `g == 0`)
#' across time periods. It returns a data frame with time periods and corresponding outcome means.
#'
#' @param y_numeric A numeric vector representing the outcome variable.
#' @param g A vector indicating group membership (e.g., treatment cohort). Groups with `g == 0` are treated as controls and excluded.
#' @param type_name_treated A character string labeling the resulting outcome series; defaults to `"Treated"`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{t}{Time period.}
#'   \item{y}{Weighted mean of the outcome variable for treated units.}
#'   \item{type}{Label identifying the outcome series.}
#' }
#'
#' @details
#' This function uses uniform weights (i.e., all weights equal to 1) and excludes observations where `g == 0`
#' (typically representing untreated or never-treated units).
#'
#' @seealso [get_weighted_means_by_gjt()]
#'
#' @export
get_df_outcome_treated <- function(y_numeric, g,
                                   type_name_treated = "Treated") {
  get_weighted_means_by_gjt(y_numeric, g, rep(1, length(g))) %>%
    rename(y = weighted_mean) %>%
    filter(g != 0) %>%
    ungroup() %>%
    select(t, y) %>%
    mutate(type = type_name_treated)
}

#' Construct Counterfactual Outcomes Using a Transition Model
#'
#' Computes a counterfactual outcome trajectory based on a transition model by subtracting estimated ATT
#' values from the observed outcomes for treated units. The output is aligned with the treated period.
#'
#' @param transition_model A fitted transition model object containing elements such as `Ps_control`,
#'   `treated_period`, and a compatible input for `summarize_transition_model()`.
#' @param y_numeric A numeric vector of observed outcomes.
#' @param g A vector of treatment group indicators (e.g., cohort indicators).
#' @param display_J Logical; if `TRUE`, includes the `J` value in the `type` label. Default is `FALSE`.
#' @param type_name_counterfactual_ours A character string for labeling the counterfactual series. Default is `"Counterfactual Untreated"`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{t}{Time period.}
#'   \item{y}{Estimated untreated (counterfactual) outcome.}
#'   \item{type}{Label identifying the counterfactual series.}
#' }
#'
#' @seealso [summarize_transition_model()], [get_df_outcome_treated()]
#'
#' @export
get_df_counterfactual_transition_model <- function(transition_model, y_numeric, g, display_J = FALSE,
                                                   type_name_counterfactual_ours = "Our Method") {
  # Extract ATT values for both models
  J <- length(transition_model$priors)
  T_max <- length(transition_model$Ps_control[[1]]) + 1
  treated_period <- transition_model$treated_period
  model_summary <- summarize_transition_model(transition_model)
  att_dbt <- c(0, sapply(list(model_summary), "[[", "ATT"))

  # Compute df_outcome_treated
  df_outcome_treated <- get_df_outcome_treated(y_numeric, g)

  # Time periods for counterfactuals
  ts <- (T_max - length(att_dbt) + 1):T_max

  # Construct counterfactual dataset
  type_name_counterfactual_ours <- paste0(type_name_counterfactual_ours,
                                          ifelse(display_J, paste0(", J = ", J), ""))
  data.frame(att = att_dbt,
             t = ts,
             type = rep(type_name_counterfactual_ours, each = length(ts))) %>%
    left_join(df_outcome_treated, by = c("t" = "t"), suffix = c("", ".y")) %>%
    mutate(untreated_counterfactual = y - att) %>%
    mutate(y = untreated_counterfactual) %>%
    select(t, y, type) %>%
    filter(t >= treated_period - 1)
}

#' Construct Counterfactual Outcomes Using a Difference-in-Differences Model
#'
#' Computes a counterfactual outcome trajectory based on a DiD model by subtracting estimated ATT values
#' from the observed treated outcomes. The result estimates untreated outcomes over time.
#'
#' @param did_model A list containing ATT estimates (`att`) and treatment group indicators (`group`).
#' @param y_numeric A numeric vector of observed outcomes.
#' @param g A vector of treatment group indicators.
#' @param type_name_counterfactual_did A character string for labeling the DiD counterfactual series. Default is `"Counterfactual Untreated (DiD)"`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{t}{Time period.}
#'   \item{y}{Estimated untreated (counterfactual) outcome.}
#'   \item{type}{Label identifying the counterfactual series.}
#' }
#'
#' @seealso [get_df_outcome_treated()]
#'
#' @export
get_df_counterfactual_did_model <- function(did_model, y_numeric, g,
  type_name_counterfactual_did = "DiD") {
  # Extract ATT values for both models
  T_max <- length(did_model$att) + 1
  treated_period <- max(did_model$group)
  att_did <- c(0, did_model$att[(treated_period - 1):(T_max - 1)])

  # Compute df_outcome_treated
  df_outcome_treated <- get_df_outcome_treated(y_numeric, g)

  # Time periods for counterfactuals
  ts <- (T_max - length(att_did) + 1):T_max

  # Construct counterfactual dataset
  data.frame(att = att_did,
             t = ts,
             type = rep(type_name_counterfactual_did, each = length(ts))) %>%
    left_join(df_outcome_treated, by = c("t" = "t"), suffix = c("", ".y")) %>%
    mutate(untreated_counterfactual = y - att) %>%
    mutate(y = untreated_counterfactual) %>%
    select(t, y, type) %>%
    filter(t >= treated_period - 1)
}

#' Get Data Frame of Counterfactuals for More Lag Model
#'
#' This function retrieves the counterfactuals for a more lagged model.
#' @param did_conditional_lag_model The more lagged model object.
#' @param y_numeric The numeric outcome data.
#' @param g The treatment assignment vector.
#' @return A data frame containing the counterfactuals for the more lagged model.
get_df_counterfactual_did_conditional_lag_model <- function(did_conditional_lag_model, y_numeric, g) {
  lag <- did_conditional_lag_model$lag
  type_string <- paste("Our Method, J = 1 (k =", lag, "Lags)")
  if (lag == 1) {
    type_string <- "Our Method, J = 1"
  }

  # Get the counterfactuals for the did_conditional_lag_model
  get_df_counterfactual_did_model(did_conditional_lag_model, y_numeric, g, type_string)
}

#' Get Data Frame of Counterfactuals for More Lag Models
#'
#' This function retrieves the counterfactuals for a list of more lagged models.
#' @param did_conditional_lag_models A list of more lagged model objects.
#' @param y_numeric The numeric outcome data.
#' @param g The treatment assignment vector.
#' @return A data frame containing the counterfactuals for all more lagged models.
get_df_counterfactual_did_conditional_lag_models <- function(did_conditional_lag_models, y_numeric, g) {
  # Initialize an empty list to store data frames
  df_list <- lapply(did_conditional_lag_models, get_df_counterfactual_did_conditional_lag_model,
    y_numeric = y_numeric, g = g)

  # Combine all data frames into one
  combined_df <- dplyr::bind_rows(df_list)
  return(combined_df)
}

#' @title Create ATT Summary Data Frame from a DiD Model
#' @description Extracts ATT (average treatment effect on the treated) estimates
#'   and 95 percent confidence intervals from a \code{did} model object. ATT
#'   estimates are taken from \code{did_model$att} and standard errors from
#'   \code{did_model$se}. Confidence intervals use the normal approximation.
#' @param did_model A model object estimated using the \code{did} package,
#'   containing vectors \code{att}, \code{se}, and \code{group}.
#' @param type_name A character string to label the output data series.
#' @return A data frame with columns \code{t} (time period), \code{estimate}
#'   (ATT point estimate), \code{ci_lower}, \code{ci_upper}, and \code{type}.
#' @seealso \code{\link[did]{att_gt}}
#' @keywords internal
get_df_att_did_model <- function(did_model, type_name = "DiD") {
  # Extract ATT values from DiD model from TWFE
  T_max <- length(did_model$att)
  treated_period <- max(did_model$group)
  att_did <- did_model$att[(treated_period):T_max]
  att_did_se <- did_model$se[(treated_period):T_max]

  # Extract ATT values from DiD model if estimated from a did package
  if (class(did_model) == "MP") {
    T_max <- length(did_model$att) + 1
    treated_period <- max(did_model$group)
    att_did <- did_model$att[(treated_period - 1):(T_max - 1)]
    att_did_se <- did_model$se[(treated_period - 1):(T_max - 1)]
  }

  # Construct data frame
  att_did_upper <- att_did + qnorm(0.975)*att_did_se
  att_did_lower <- att_did - qnorm(0.975)*att_did_se
  ts <-  1:(T_max-treated_period + 1)
  data.frame(t = ts,
             estimate = att_did,
             ci_lower = att_did_lower,
             ci_upper = att_did_upper,
             type = type_name)
}

#' Create ATT Data Frame from Conditional Lag DiD Model
#'
#' Extracts ATT estimates from a conditional lag DiD model and formats them
#' as a data frame with appropriate type labels.
#'
#' @param did_conditional_lag_model A conditional lag DiD model object
#'   containing a lag field.
#'
#' @return A data frame with ATT estimates and confidence intervals.
#' @keywords internal
get_df_att_did_conditional_lag_model <- function(did_conditional_lag_model) {
  lag <- did_conditional_lag_model$lag
  type_string <- paste("Our Method, J = 1 (k =", lag,  "Lags)")
  if (lag == 1) {
    type_string <- "Our Method, J = 1"
  }
  get_df_att_did_model(did_conditional_lag_model, type_string)
}

#' Get Data Frame of ATTs for More Lag Models
#'
#' This function retrieves the ATTs for a list of more lagged models.
#' @param did_conditional_lag_models A list of more lagged model objects.
#' @return A data frame containing the ATTs for all more lagged models.
#' @keywords internal
get_df_att_did_conditional_lag_models <- function(did_conditional_lag_models) {
  # Initialize an empty list to store data frames
  df_list <- lapply(did_conditional_lag_models, get_df_att_did_conditional_lag_model)

  # Combine all data frames into one
  combined_df <- dplyr::bind_rows(df_list)
  return(combined_df)
}


#' @title Create ATT Summary Data Frame from a Transition Model
#' @description Extracts ATT (average treatment effect on the treated) estimates
#'   and confidence intervals from a fitted transition model. Point estimates
#'   come from \code{summarize_transition_model()} and confidence intervals from
#'   bootstrap estimates if available.
#' @param transition_model A transition model object containing fields such as
#'   \code{Ps_control_empirical}, \code{treated_period}, \code{priors}, and
#'   \code{bootstrap_estimates}.
#' @param param_name Name of the parameter to extract (default \code{"ATT"}).
#' @param uniform_ci Logical; if \code{TRUE}, compute uniform (simultaneous)
#'   confidence intervals; if \code{FALSE}, compute pointwise intervals.
#' @return A data frame with columns \code{t} (time period), \code{estimate},
#'   \code{ci_lower}, \code{ci_upper}, and \code{type}.
#' @seealso \code{\link{summarize_transition_model}}
#' @keywords internal
get_df_att_transition_model <- function(transition_model, param_name = "ATT",
                                        uniform_ci = FALSE) {
  # Prepare data
  T_max <- length(transition_model$Ps_control_empirical[[1]]) + 1
  treated_period <- transition_model$treated_period
  J <- length(transition_model$priors)
  model_summary <- summarize_transition_model(transition_model)

  # Extract ATT estimates and uniform confidence intervals
  estimate <- model_summary[[param_name]]
  t <- names(estimate)
  if(is.null(names(estimate)[1])) {
    t <- 1:length(estimate)
    names(estimate) <- t
  }
  ci_upper <- NA
  ci_lower <- NA
  bootstrap_models <- transition_model$bootstrap_estimates

  df_params <- data.frame(t = t, estimate = c(estimate))

  if (!is.null(bootstrap_models)) {
    bootstrap_models_summary <- summarize_transition_models(bootstrap_models)
    bootstrap_estimates <- lapply(bootstrap_models_summary, "[[", param_name)

    df_params <- bootstrap_estimates_to_ci_df(estimate, bootstrap_estimates,
                                              uniform_ci = uniform_ci)

    has_ci <- TRUE
  }

  # Form a data frame
  J_names <- paste0("Our Method, J = ", J)

  # Attach CIs accordingly
  df_params %>%
    mutate(type = J_names,
           ci_upper = as.numeric(ci_upper),
           ci_lower = as.numeric(ci_lower)) %>%
    select(t, estimate, ci_lower, ci_upper, type)
}

