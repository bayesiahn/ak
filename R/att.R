#' Calculate the Difference in Outcomes by Transition Probabilities
#'
#' This function computes the difference in expected outcomes between treated and
#' control groups based on their respective transition probabilities and the
#' distribution of initial states.
#'
#' @param pmf A probability mass function (vector) representing the distribution
#'        of initial states.
#' @param P_treated A transition matrix for the treated group.
#' @param P_control A transition matrix for the control group.
#' @param Y_values A vector of outcome values.
#'
#' @return A numeric value representing the difference in expected outcomes
#'         between the treated and control groups.
#'
#' @examples
#' pmf <- c(0.5, 0.5)
#' P_treated <- matrix(c(0.7, 0.4, 0.3, 0.6), nrow = 2)
#' P_control <- matrix(c(0.6, 0.5, 0.4, 0.5), nrow = 2)
#' Y_values <- c(0, 1)
#' difference <- difference_by_transitions(pmf, P_treated, P_control, Y_values)
#'
#' @export
difference_by_transitions <- function(pmf, P_treated, P_control, Y_values) {
  get_expected_future_outcome(pmf, P_treated - P_control, Y_values)
}

#' Calculate Expected Future Outcome
#'
#' This function computes the expected future outcome given an initial distribution
#' of states (pmf), a future transition matrix (P_future), and outcome values (Y_values).
#'
#' @param pmf A probability mass function (vector) representing the distribution
#'        of initial states.
#' @param P_future A transition matrix representing future state transitions.
#'        This can be the difference between treated and control transition matrices
#'        to compute the difference in expected outcomes.
#' @param Y_values A vector of outcome values corresponding to each state.
#'
#' @return A numeric value representing the expected future outcome.
#'
#' @examples
#' pmf <- c(0.5, 0.5)
#' P_future <- matrix(c(0.1, -0.1, -0.1, 0.1), nrow = 2)
#' Y_values <- c(0, 1)
#' expected_outcome <- get_expected_future_outcome(pmf, P_future, Y_values)
#' print(expected_outcome)
#'
#' @export
get_expected_future_outcome <- function(pmf, P_future, Y_values) {
  get_expected_outcome(t(P_future) %*% pmf, Y_values)
}

#' Calculate Expected Outcome
#'
#' Computes the expected outcome given a probability mass function (pmf) and
#' a vector of outcome values (Y_values).
#'
#' @param pmf A probability mass function (vector) representing the distribution
#'        of initial states.
#' @param Y_values A vector of outcome values.
#'
#' @return A numeric value representing the expected outcome based on the pmf and Y_values.
#'
#' @examples
#' pmf <- c(0.5, 0.5)
#' Y_values <- c(0, 1)
#' expected_outcome <- get_expected_outcome(pmf, Y_values)
#' print(expected_outcome)
#'
#' @export
get_expected_outcome <- function(pmf, Y_values) {
  t(pmf) %*% Y_values
}

#' Expected Outcomes from One‑Step Transitions (composed recursively)
#'
#' Compute expected outcomes over time when you have a sequence of **one‑step**
#' transition matrices (e.g., period‑specific Markov transitions). For each
#' period \eqn{t}, the function composes \eqn{P_1 P_2 \cdots P_t} internally and
#' applies it to the initial distribution to obtain the distribution at \eqn{t}.
#'
#' @param pmf_initial Numeric probability vector of initial state distribution
#'   (length = number of outcome states). Must sum to 1.
#' @param Ps List of **one‑step** transition matrices for periods 1..T.
#'   Each \code{Ps[[t]]} is a square stochastic matrix whose dimensions match
#'   \code{length(pmf_initial)} and \code{length(Y_values)}.
#' @param Y_values Numeric vector of outcome values (same order as state space).
#'
#' @return Numeric vector of length \eqn{T} with expected outcomes at each period.
#'
#' @details
#' Use this function when your inputs are per‑period transitions and you want the
#' function to handle the cumulative composition \eqn{P^{(t)} = \prod_{s=1}^t P_s}.
#' This is the safest choice when transitions vary by period.
#'
#' @note
#' This function performs matrix products up to \eqn{t} for each period, which can
#' be more computationally expensive than using pre‑computed \eqn{t}-step matrices.
#'
#' @seealso
#' \code{\link{calculate_expected_outcomes}} for the version that expects
#' pre‑computed \eqn{t}-step transition matrices; \code{\link{multiply_matrices_in_list}}
#' which is used internally for composition.
#'
#' @examples
#' \dontrun{
#' # One-step transitions (vary by period)
#' pmf0 <- c(0.6, 0.4)
#' P1 <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE)
#' P2 <- matrix(c(0.85,0.15,0.25,0.75), 2, byrow = TRUE)
#' Ys <- c(0, 1)
#' calculate_expected_outcomes_recursively(pmf0, list(P1, P2), Ys)
#' }
#' @export
calculate_expected_outcomes_recursively <- function(pmf_initial, Ps, Y_values) {
  T_max <- length(Ps)
  expected_outcomes <- rep(0, T_max)

  # Use cumulative products to avoid recomputation: O(T) instead of O(T²)
  P_cumulative <- NULL
  for (t in 1:T_max) {
    # Update cumulative product incrementally
    if (t == 1) {
      P_cumulative <- Ps[[1]]
    } else {
      P_cumulative <- P_cumulative %*% Ps[[t]]
    }
    # Expected outcome at t
    expected_outcomes[t] <- sum(pmf_initial %*% P_cumulative %*% Y_values)
  }

  return(expected_outcomes)
}

#' Expected Outcomes from Pre‑computed t‑Step Transitions (no recursion)
#'
#' Compute expected outcomes over time when you already have the **cumulative**
#' \eqn{t}-step transition matrix for each period \eqn{t}. For each \eqn{t}, the
#' function applies \eqn{P^{(t)}} directly to the initial distribution.
#'
#' @param pmf_initial Numeric probability vector of initial state distribution
#'   (length = number of outcome states). Must sum to 1.
#' @param Ps List of **cumulative** \eqn{t}-step transition matrices for periods 1..T.
#'   Each \code{Ps[[t]]} must equal \eqn{P^{(t)} = \prod_{s=1}^t P_s} and be a
#'   square stochastic matrix consistent with \code{pmf_initial} and \code{Y_values}.
#' @param Y_values Numeric vector of outcome values (same order as state space).
#'
#' @return Numeric vector of length \eqn{T} with expected outcomes at each period.
#'
#' @details
#' Choose this function if you’ve pre‑computed \eqn{P^{(t)}} (e.g., for speed or
#' numerical control) and want to avoid repeated in‑function matrix multiplications.
#' If you only have one‑step transitions, use
#' \code{\link{calculate_expected_outcomes_recursively}} instead.
#'
#' @note
#' Passing one‑step transitions here will yield incorrect results. Ensure that
#' \code{Ps[[t]]} is the **t‑step** transition matrix from the initial period to \eqn{t}.
#'
#' @seealso
#' \code{\link{calculate_expected_outcomes_recursively}} for composing one‑step
#' transitions on the fly.
#'
#' @examples
#' \dontrun{
#' # Pre-computed t-step transitions
#' pmf0 <- c(0.6, 0.4); Ys <- c(0, 1)
#' P1 <- matrix(c(0.9,0.1,0.2,0.8), 2, byrow = TRUE)     # 1-step
#' P2 <- P1 %*% matrix(c(0.85,0.15,0.25,0.75), 2, TRUE)  # 2-step (pre-computed)
#' calculate_expected_outcomes(pmf0, list(P1, P2), Ys)
#' }
#' @export
calculate_expected_outcomes <- function(pmf_initial, Ps, Y_values) {
  T_max <- length(Ps)
  expected_outcomes <- rep(0, T_max)

  for (t in 1:T_max) {
    # Ps[[t]] must be the t-step transition matrix
    expected_outcomes[t] <- sum(pmf_initial %*% Ps[[t]] %*% Y_values)
  }

  return(expected_outcomes)
}

#' LTATT from One-Step Markov Transitions (recursive composition)
#'
#' Compute latent-type ATTs (LTATTs) over post-treatment periods by differencing
#' expected outcomes for treated and control transition regimes. This version expects
#' **period-specific one-step** transition matrices for each group and composes them
#' internally up to period \eqn{t} (i.e., uses \eqn{P_1 P_2 \cdots P_t}).
#'
#' @param pmf_treated_just_before_treatment Numeric probability vector (sums to 1)
#'   for the treated group's state distribution at the last pre-treatment period.
#' @param Ps_treated List of **one-step** transition matrices for the treated regime,
#'   ordered by post-treatment period (same square dimension as the state space).
#' @param Ps_control List of **one-step** transition matrices for the control regime,
#'   ordered by post-treatment period; must have the same length as \code{Ps_treated}.
#' @param Y_values Numeric vector of outcome values (aligned with the state order).
#'
#' @return Numeric vector of LTATTs for each post-treatment period.
#'
#' @details
#' Use this function when your inputs are **per-period one-step** Markov transitions and
#' you want the function to handle the recursive composition into \eqn{t}-step transitions.
#' If you already have \eqn{t}-step (cumulative) transitions, use \code{\link{compute_LTATTs}} instead.
#'
#' @note Passing cumulative \eqn{t}-step matrices here will effectively “double count”
#' transitions and yield incorrect LTATTs. Ensure \code{Ps_*} are one-step matrices.
#'
#' @seealso \code{\link{compute_LTATTs}},
#'   \code{\link{calculate_expected_outcomes_recursively}},
#'   \code{\link{estimate_post_em_algorithm_parameters}}.
#'
#' @examples
#' \dontrun{
#' pmf0 <- c(0.7, 0.3); Ys <- c(0, 1)
#' P1_tr <- matrix(c(0.9,0.1, 0.2,0.8), 2, byrow=TRUE)
#' P2_tr <- matrix(c(0.85,0.15, 0.25,0.75), 2, byrow=TRUE)
#' P1_co <- matrix(c(0.92,0.08, 0.15,0.85), 2, byrow=TRUE)
#' P2_co <- matrix(c(0.90,0.10, 0.20,0.80), 2, byrow=TRUE)
#' compute_LTATTs_under_markovian(pmf0, list(P1_tr,P2_tr), list(P1_co,P2_co), Ys)
#' }
#' @export
compute_LTATTs_under_markovian <- function(pmf_treated_just_before_treatment, Ps_treated, Ps_control, Y_values) {
  # Check if the length of treated and control transition matrices lists are the same
  if (length(Ps_treated) != length(Ps_control)) {
    stop("The number of transition matrices for the treated group and the control group must be the same.")
  }

  # Calculate expected outcomes for treated and control groups
  expected_treated <- calculate_expected_outcomes_recursively(pmf_treated_just_before_treatment, Ps_treated, Y_values)
  expected_control <- calculate_expected_outcomes_recursively(pmf_treated_just_before_treatment, Ps_control, Y_values)

  # Compute LTATTs as the difference between treated and control outcomes
  LTATTs <- expected_treated - expected_control

  return(LTATTs)
}

#' LTATT from Pre-computed \eqn{t}-Step Counterfactual Transitions (no recursion)
#'
#' Compute latent-type ATTs (LTATTs) over post-treatment periods using
#' pre-computed **\eqn{t}-step (cumulative)** transition matrices that encode the
#' counterfactual evolution “as if” post-treatment dynamics followed pre-treatment
#' patterns (e.g., outputs of \code{estimate_post_em_algorithm_parameters()}).
#'
#' @param pmf_treated_just_before_treatment Numeric probability vector (sums to 1)
#'   for the treated group's state distribution at the last pre-treatment period.
#' @param Ps_treated_from_pre List of **cumulative \eqn{t}-step** transition matrices
#'   for the treated group’s counterfactual path implied by pre-treatment dynamics.
#'   Each \code{Ps_treated_from_pre[[t]]} should equal \eqn{P^{(t)} = \prod_{s=1}^t P_s}.
#' @param Ps_control_from_pre List of **cumulative \eqn{t}-step** transition matrices
#'   for the control group’s path implied by pre-treatment dynamics; must have the
#'   same length as \code{Ps_treated_from_pre}.
#' @param Y_values Numeric vector of outcome values (aligned with the state order).
#'
#' @return Numeric vector of LTATTs for each post-treatment period.
#'
#' @details
#' This function implements the transition-independence counterfactual used to identify
#' LTATT/ATT in event-study designs with discrete outcomes: expected untreated outcomes
#' for treated units are constructed by applying (latent-type-specific) pre-treatment
#' transition dynamics forward. See Ahn & Kasahara (2025) for the framework and estimands.
#'
#' @note Passing one-step per-period matrices here will yield incorrect results. Ensure
#' \code{Ps_*_from_pre[[t]]} are **\eqn{t}-step** transitions. If you only have one-step
#' transitions, use \code{\link{compute_LTATTs_under_markovian}} instead.
#'
#' @seealso \code{\link{compute_LTATTs_under_markovian}},
#'   \code{\link{calculate_expected_outcomes}},
#'   \code{\link{estimate_post_em_algorithm_parameters}}.
#'
#' @examples
#' \dontrun{
#' pmf0 <- c(0.7, 0.3); Ys <- c(0, 1)
#' # Suppose P1 and P2_pre are 1-step and 2-step transitions implied from pre-treatment
#' P1 <- matrix(c(0.9,0.1, 0.2,0.8), 2, byrow=TRUE)
#' P2_pre <- P1 %*% matrix(c(0.85,0.15, 0.25,0.75), 2, byrow=TRUE)  # 2-step cumulative
#' compute_LTATTs(pmf0, list(P1, P2_pre), list(P1, P2_pre), Ys)
#' }
#' @export
compute_LTATTs <- function(pmf_treated_just_before_treatment, Ps_treated_from_pre, Ps_control_from_pre, Y_values) {
  # Check if the length of treated and control transition matrices lists are the same
  if (length(Ps_treated_from_pre) != length(Ps_control_from_pre)) {
    stop("The number of transition matrices for the treated group and the control group must be the same.")
  }

  # Calculate expected outcomes for treated and control groups
  expected_treated <- calculate_expected_outcomes(pmf_treated_just_before_treatment, Ps_treated_from_pre, Y_values)
  expected_control <- calculate_expected_outcomes(pmf_treated_just_before_treatment, Ps_control_from_pre, Y_values)

  # Compute LTATTs as the difference between treated and control outcomes
  LTATTs <- expected_treated - expected_control

  return(LTATTs)
}
