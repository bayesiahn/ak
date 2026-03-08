#' Difference-in-Differences with Conditional Parallel Trends
#'
#' This function calculates the difference-in-differences (DiD) estimates using conditional parallel trends.
#' It uses the `get_conditional_state_matrix` function to compute lagged values of the outcome variable `y`.
#'
#' @param y A numeric matrix of dimensions N by T, representing the outcome variable.
#' @param g A numeric vector of length N, where each element is zero if not treated and is equal to `T_treated` if treated in that particular period.
#' @param lags An integer specifying the number of lags for the conditional states. Default is 1.
#'
#' @return A numeric value representing the difference-in-differences estimate.
#'
#' @examples
#' y <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3, byrow = TRUE)
#' g <- c(0, 0, 2)
#' did_conditional_parallel_trends(y, g, lags = 1)
#'
#' @export
did_conditional_parallel_trends <- function(y, g, lags = 1) {
  T_treated <- unique(g[g > 0])
  if (length(T_treated) != 1) {
    stop("There should be exactly one common treatment period `T_treated`.")
  }
  T_treated <- T_treated[1]

  N <- nrow(y)
  T <- ncol(y)

  # Get conditional state matrix
  conditional_matrix <- get_conditional_state_matrix(y, lags)

  # Treated and control groups
  treated_indices <- which(g == T_treated)
  control_indices <- which(g == 0)

  conditional_matrix_treated <- matrix(conditional_matrix[treated_indices, ], ncol = T)
  conditional_matrix_control <- matrix(conditional_matrix[control_indices, ], ncol = T)
  y_treated <- matrix(y[treated_indices, ], ncol = T)
  y_control <- matrix(y[control_indices, ], ncol = T)

  # Pre-treatment period
  pre_treatment_period <- T_treated - 1
  post_treatment_period <- T_treated

  if (pre_treatment_period < 1 || pre_treatment_period <= lags) {
    stop("Pre-treatment period is too short to compute the lags.")
  }

  # Averages for each unique state and get their distribution in pre_treatment_period
  treated_pre_treatment_distribution <- table(conditional_matrix[treated_indices, pre_treatment_period]) / length(treated_indices)
  control_pre_treatment_distribution <- table(conditional_matrix[control_indices, pre_treatment_period]) / length(control_indices)
  treated_states <- names(treated_pre_treatment_distribution)

  # Get average for each unique state
  CATT <- matrix(NA, nrow = length(treated_states), ncol = T)

  for (k in 1:length(treated_states)) {
    treated_state <- treated_states[k]
    treated_unit_indices_k <- which(conditional_matrix_treated[, pre_treatment_period] == treated_state)
    control_unit_indices_k <- which(conditional_matrix_control[, pre_treatment_period] == treated_state)
    for (t in post_treatment_period:T) {
      CATT[k, t] <- mean(y_treated[treated_unit_indices_k, t]) - mean(y_control[control_unit_indices_k, t])
    }
  }

  ATT <- apply(CATT, 2, function(col_t) sum(col_t * treated_pre_treatment_distribution))
  ATT <- ATT[post_treatment_period:T]
  CATTs <- lapply(1:nrow(CATT), function(i) CATT[i, post_treatment_period:T])

  # Save as transition_model_estimate
  transition_model_estimate <- list(ATT = ATT, LTATTs = list(matrix(ATT)), CATT = CATT, treated_states = treated_states,
                              treated_pre_treatment_distribution = treated_pre_treatment_distribution,
                              control_pre_treatment_distribution = control_pre_treatment_distribution,
                              transition_model = list(priors = 1, priors_treated = 1, treated_period = T_treated),
                              weights = matrix(rep(1, N), nrow = N),
                              posteriors = matrix(rep(1, N), nrow = N))
  class(transition_model_estimate) <- "TransitionModel"

  return(transition_model_estimate)
}
