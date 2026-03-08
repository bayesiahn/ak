#' Split State Sequences by Treatment Status and Period
#'
#' This function processes a matrix of observed state sequences (`y_indices_matrix`)
#' and splits it according to treatment status (`g`) and treatment timing. It also computes
#' indicators used in the E- and M-steps of the EM algorithm, including initial outcome
#' and transition indicators.
#'
#' @param y_indices_matrix A matrix where each row corresponds to a unit and each column
#'        represents a time period, with values indicating observed outcomes or states.
#' @param g A numeric vector indicating treatment status for each unit. Units with \code{g[i] > 0}
#'        are treated in period \code{g[i]}.
#' @param e_step_control_pre_only Logical; if \code{TRUE}, the E-step for control units
#'   uses only pre-treatment observations (i.e., periods strictly before the treatment
#'   period). If \code{FALSE}, all available control-unit observations are used.
#'   Default is \code{FALSE}.
#' @param e_step_treated_pre_only Logical; if \code{TRUE}, the E-step for treated units
#'   uses only pre-treatment observations. If \code{FALSE}, all available treated-unit
#'   observations are used. Default is \code{FALSE}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{y_indices_matrix}{Original outcome/state matrix.}
#'   \item{g}{Treatment status vector.}
#'   \item{y_indices_pretreatment}{Combined pre-treatment sequences for all units.}
#'   \item{y_indices_post_treatment}{Combined post-treatment sequences for all units.}
#'   \item{y_indices_control}{State sequences for control units.}
#'   \item{y_indices_treated}{State sequences for treated units.}
#'   \item{y_indices_control_pretreatment}{Pre-treatment sequences for control units.}
#'   \item{y_indices_treated_pretreatment}{Pre-treatment sequences for treated units.}
#'   \item{y_indices_control_post_treatment}{Post-treatment sequences for control units.}
#'   \item{y_indices_treated_post_treatment}{Post-treatment sequences for treated units.}
#'   \item{y_indices_for_e_step_control}{State sequences for control units used in E-step.}
#'   \item{y_indices_for_e_step_treated}{State sequences for treated units used in E-step.}
#'   \item{initial_outcome_indicators}{Initial outcome indicators for all units.}
#'   \item{initial_outcome_indicators_control}{Subset of indicators for control units.}
#'   \item{initial_outcome_indicators_treated}{Subset of indicators for treated units.}
#'   \item{initial_outcome_indicators_joint}{Initial outcome indicators for joint distribution for (Y1, D); first half columns represent ones for D = 0.}
#'   \item{transition_indicators_pretreatment}{Transition indicators for pre-treatment periods.}
#'   \item{transition_indicators_control_post}{Transition indicators for control units in post-treatment periods.}
#'   \item{transition_indicators_treated_post}{Transition indicators for treated units in post-treatment periods.}
#'   \item{transition_indicators_control_post_baseline}{Transition indicators for control units in post-treatment periods, fixed at period just before treatment.}
#'   \item{transition_indicators_treated_post_baseline}{Transition indicators for treated units in post-treatment periods, fixed at period just before treatment.}
#'   \item{is_treated}{Logical vector indicating treated units.}
#'   \item{is_control}{Logical vector indicating control units.}
#'   \item{n}{Total number of units.}
#'   \item{n_control}{Number of control units.}
#'   \item{n_treated}{Number of treated units.}
#'   \item{prob_treated}{Proportion of treated units.}
#'   \item{treated_period}{Period when treatment starts for treated units.}
#' }
#'
#' @examples
#' \dontrun{
#' y_matrix <- matrix(sample(1:3, 20, replace = TRUE), nrow = 5)
#' g <- c(0, 1, 0, 1, 0)
#' data_for_est <- split_data_for_estimation(y_matrix, g)
#' }
split_data_for_estimation <- function(y_indices_matrix, g,
                                      e_step_control_pre_only = FALSE,
                                      e_step_treated_pre_only = FALSE) {
  is_treated <- g > 0
  is_control <- g <= 0
  n_control = sum(is_control)
  n_treated = sum(is_treated)
  prob_treated <- n_treated / (n_control + n_treated)

  if (!n_treated > 0) {
    stop("There are no treated units available.")
  }
  if (ncol(y_indices_matrix) < 2) {
    stop("y_indices must have at least two columns.")
  }

  treated_period <- ifelse(n_treated > 0, min(g[is_treated]), 0)
  T_max <- ncol(y_indices_matrix)
  just_before_treatment_periods <- (treated_period-1):(treated_period-1)
  y_indices_just_before_treatment_periods <- matrix(y_indices_matrix[,just_before_treatment_periods],
                                    ncol = length(just_before_treatment_periods))

  K <- length(unique(c(y_indices_matrix)))
  supports <- 1:K

  # Indicators
  initial_outcome_indicators <- y_indices_to_initial_outcome_indicators(y_indices_matrix, supports)
  just_before_treatment_outcome_indicators <- as.matrix(y_indices_to_initial_outcome_indicators(y_indices_just_before_treatment_periods, supports))
  initial_outcome_indicators_control <- matrix(initial_outcome_indicators[is_control,], nrow = n_control)
  initial_outcome_indicators_treated <- matrix(initial_outcome_indicators[is_treated,], nrow = n_treated)
  initial_outcome_indicators_joint <- Matrix::bdiag(initial_outcome_indicators_control, initial_outcome_indicators_treated)

  indicators_treated_just_before_treatment <- matrix(just_before_treatment_outcome_indicators[is_treated,],
      nrow = n_treated)
  indicators_control_just_before_treatment <- matrix(just_before_treatment_outcome_indicators[is_control,],
      nrow = n_control)

  # Split the outcome/state matrix by treatment status
  y_indices_control <- matrix(y_indices_matrix[is_control, ], nrow = n_control)
  y_indices_treated <- matrix(y_indices_matrix[is_treated, ], nrow = n_treated)

  y_indices_control_pretreatment <- NA
  y_indices_treated_pretreatment <- NA
  if (treated_period > 1) {
    y_indices_control_pretreatment <- matrix(y_indices_control[, 1:(treated_period - 1)],
                                             nrow = n_control)
    y_indices_treated_pretreatment <- matrix(y_indices_treated[, 1:(treated_period - 1)],
                                             nrow = n_treated)
  }
  y_indices_control_post_treatment <- matrix(y_indices_control[, (treated_period-1):ncol(y_indices_control)],
                                              nrow = n_control)
  y_indices_treated_post_treatment <- matrix(y_indices_treated[, (treated_period-1):ncol(y_indices_treated)],
                                            nrow = n_treated)


  # Combine pre- and post-treatment data
  y_indices_pretreatment <- rbind(y_indices_control_pretreatment, y_indices_treated_pretreatment)
  y_indices_post_treatment <- rbind(y_indices_control_post_treatment, y_indices_treated_post_treatment)

  # Prepare data for E-step
  y_indices_for_e_step_control <- y_indices_control
  y_indices_for_e_step_treated <- y_indices_treated

  # Use only pre-treatment data for control and treated units
  if (e_step_control_pre_only) {
    y_indices_for_e_step_control <- y_indices_control_pretreatment
  }
  if (e_step_treated_pre_only) {
    y_indices_for_e_step_treated <- y_indices_treated_pretreatment
  }

  # Transition indicators
  transition_indicators <- transitions_to_transition_indicators(y_indices_matrix, supports)
  transition_indicators_pretreatment <- transitions_to_transition_indicators(y_indices_pretreatment, supports)
  transition_indicators_control <- transitions_to_transition_indicators(y_indices_control, supports)
  transition_indicators_treated <- transitions_to_transition_indicators(y_indices_treated, supports)

  # Transition indicators for post-treatment periods
  transition_indicators_control_post <- transitions_to_transition_indicators(y_indices_control_post_treatment, supports)
  transition_indicators_treated_post <- transitions_to_transition_indicators(y_indices_treated_post_treatment, supports)

  # Transition indicators for post-treatment periods, but fixing period to be fixed at 1
  transition_indicators_control_post_baseline <- transitions_to_transition_indicators(y_indices_control_post_treatment, supports,
                                                                                      current_period_to_be_fixed = 1)
  transition_indicators_treated_post_baseline <- transitions_to_transition_indicators(y_indices_treated_post_treatment, supports,
                                                                                      current_period_to_be_fixed = 1)

  # Transition indicators, but fixing period to be fixed at 1
  transition_indicators_control_baseline <- transitions_to_transition_indicators(y_indices_control, supports,
                                                                                 current_period_to_be_fixed = 1)
  transition_indicators_treated_baseline <- transitions_to_transition_indicators(y_indices_treated, supports,
                                                                                 current_period_to_be_fixed = 1)

  list(
    y_indices_matrix = y_indices_matrix,
    g = g,
    y_indices_pretreatment = y_indices_pretreatment,
    y_indices_post_treatment = y_indices_post_treatment,
    y_indices_control = y_indices_control,
    y_indices_treated = y_indices_treated,
    y_indices_control_pretreatment = y_indices_control_pretreatment,
    y_indices_treated_pretreatment = y_indices_treated_pretreatment,
    y_indices_control_post_treatment = y_indices_control_post_treatment,
    y_indices_treated_post_treatment = y_indices_treated_post_treatment,
    y_indices_for_e_step_control = y_indices_for_e_step_control,
    y_indices_for_e_step_treated = y_indices_for_e_step_treated,
    initial_outcome_indicators = initial_outcome_indicators,
    initial_outcome_indicators_control = initial_outcome_indicators_control,
    initial_outcome_indicators_treated = initial_outcome_indicators_treated,
    initial_outcome_indicators_joint = initial_outcome_indicators_joint,
    indicators_treated_just_before_treatment = indicators_treated_just_before_treatment,
    indicators_control_just_before_treatment = indicators_control_just_before_treatment,
    transition_indicators = transition_indicators,
    transition_indicators_pretreatment = transition_indicators_pretreatment,
    transition_indicators_control = transition_indicators_control,
    transition_indicators_treated = transition_indicators_treated,
    transition_indicators_control_post = transition_indicators_control_post,
    transition_indicators_treated_post = transition_indicators_treated_post,
    transition_indicators_control_baseline = transition_indicators_control_baseline,
    transition_indicators_treated_baseline = transition_indicators_treated_baseline,
    transition_indicators_control_post_baseline = transition_indicators_control_post_baseline,
    transition_indicators_treated_post_baseline = transition_indicators_treated_post_baseline,
    is_treated = is_treated,
    is_control = is_control,
    n = length(g),
    n_control = n_control,
    n_treated = n_treated,
    prob_treated = prob_treated,
    treated_period = treated_period
  )
}
