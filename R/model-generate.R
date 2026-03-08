#' @title Generate transition model
#' @description Generate a transition model with parameters
#'
#' @param outcomes A vector of outcomes
#' @param J Number of transition patterns
#' @param T_max Maximum number of time periods
#' @param treated_period The period in which the treatment is applied; zero or
#' negative value if never treated, which will make the model not include
#' Ps_treated
#' @export
#' @return A list containing the parameters for the transition model with J transition patterns:
#' - \code{Y}: An outcome space object
#' - \code{priors}: A J-vector of prior probabilities for each transition pattern
#' - \code{prob_treated}: A scalar representing the probability of being treated
#' - \code{priors_treated}: A J-vector of prior probabilities of belonging to each transition pattern for treated units
#' - \code{treated_period}: The period in which the treatment is applied
#' - \code{Ps_control}: A J-length list of (T_max-1)-length list of transition
#'   matrices for each transition pattern and time period for control units
#' - \code{Ps_treated}: A J-length list of (T_max-1)-length list of transition
#'   if treated_period > 0, otherwise NULL
#' - \code{Ps_control_empirical}: A J-length list of (T_max-1)-length list of transition
#'   matrices for each transition pattern and time period for control units,
#'   where the transition matrices are generated from the empirical distribution
#'   of the data
#'   \code{Ps_treated_empirical}: A J-length list of (T_max-1)-length list of transition
#'   matrices for each transition pattern and time period for treated units,
#'   where the transition matrices are generated from the empirical distribution
#'   of the data
#' - \code{pmfs_initial_control}: A J-length list of pmfs for initial outcomes for each transition pattern for control units
#' - \code{pmfs_initial_treated}: A J-length list of pmfs for initial outcomes for each transition pattern for treated units
#' - \code{pmfs_pre_treatment_control}: A J-length list of pmfs for outcomes at (treated_period-1)th period for each transition pattern for control units
#' - \code{pmfs_pre_treatment_treated}: A J-length list of pmfs for outcomes at (treated_period-1)th period for each transition pattern for treated units
#' - \code{p_y1d}: A J-length list of joint pmfs for (Y_1, D) represented by K by 2 matrices (1st column is for D = 0, 2nd column is for D = 1)
#' -
generate_transition_model <- function(outcomes = c(0,1), J = 2, T_max = 10,
                                      treated_period = 10) {
  # Generate outcome space
  Y <- generate_Y_t(outcomes)

  # Generate transition matrices
  treated_period <- min(treated_period, T_max)
  Ps_control <- generate_Ps(Y, J, T_max)

  # Generate initial outcome distributions
  pmfs_initial_control <- generate_initial_outcome_dists(Y, J)
  pmfs_initial_treated <- generate_initial_outcome_dists(Y, J)

  # Generate transition matrices for treated
  Ps_treated <- NULL
  if (treated_period > 0) {
    treated_period <- min(treated_period, T_max)
    Ps_treated <- generate_Ps(Y, J, (T_max-treated_period+2),
                              early_null_periods = treated_period-2)

    # Fill the first (treated_period-1) periods with the control transition matrices
    if (treated_period > 2) {
      for (j in 1:J) {
        Ps_treated[[j]][1:(treated_period-2)] <- Ps_control[[j]][1:(treated_period-2)]
      }
    }
  }


  # Generate priors
  priors <- pmax(0.05, generate_pmf(1:J))
  priors <- priors / sum(priors)

  # Generate priors for treatment probability for each latent group
  p_d_given_z <- round_to_decimal_point(runif(J, min = 0.05, max = 0.95))
  prob_treated <- sum(p_d_given_z * priors)
  priors_treated <- p_d_given_z * priors / prob_treated

  # Compute pmfs_pre_treatment_control, pmfs_pre_treatment_treated
  pmfs_pre_treatment_control <- pmfs_initial_control
  pmfs_pre_treatment_treated <- pmfs_initial_treated
  if (treated_period > 2) {
    for (j in 1:J) {
      pmfs_pre_treatment_control[[j]] <- apply_transition_matrices(pmfs_initial_control[[j]],
                                                                   Ps_control[[j]][1:(treated_period-2)])
      pmfs_pre_treatment_treated[[j]] <- apply_transition_matrices(pmfs_initial_treated[[j]],
                                                                   Ps_treated[[j]][1:(treated_period-2)])
    }
  }

  # Compute p_y1d
  p_y1d <- list()
  for (j in 1:J) {
    p_y1d[[j]] <- cbind(pmfs_initial_control[[j]] * (1-p_d_given_z[j]),
                        pmfs_initial_treated[[j]] * p_d_given_z[j])
  }

  # Return list
  model <- list(
    Y = Y,
    priors = priors,
    prob_treated = prob_treated,
    priors_treated = priors_treated,
    treated_period = treated_period,
    Ps_control = Ps_control,
    Ps_treated = Ps_treated,
    Ps_control_empirical = Ps_control,
    Ps_treated_empirical = Ps_treated,
    pmfs_initial_control = pmfs_initial_control,
    pmfs_initial_treated = pmfs_initial_treated,
    pmfs_pre_treatment_control = pmfs_pre_treatment_control,
    pmfs_pre_treatment_treated = pmfs_pre_treatment_treated,
    p_y1d = p_y1d
  )

  class(model) <- "TransitionModel"
  reorder_transition_model(model)
}


#' @title Generate initial outcome distributions
#' @description Generate initial outcome distributions for each transition pattern
#'
#' @param Y_t An outcome space object
#' @param J Number of transition patterns
#' matrix and initial outcome dists should be bounded below by some positive number.
#' Default value is FALSE; when TRUE, every element of transition matrices and
#' initial outcome dists are guaranteed to be strictly positive.
#'
#' @return A J-length list of pmfs for each transition pattern
generate_initial_outcome_dists <- function(Y_t, J) {
  lapply(rep(list(Y_t$names), J), generate_pmf)
}

#' @title Generate transition matrices
#' @description Generate transition matrices for each transition pattern
#'
#' @param Y_t An outcome space object
#' @param J Number of transition patterns
#' @param T_max Maximum number of time periods in which transition matrices exist
#' @param early_null_periods Number of null periods before transition matrix appears; default = 0
#' matrix and initial outcome dists should be bounded below by some positive number.
#' Default value is FALSE; when TRUE, every element of transition matrices and
#' initial outcome dists are guaranteed to be strictly positive.
#'
#' @return A J-length list of (T_max+early_null_periods-1)-length list of transition matrices for each transition pattern and time period
generate_Ps <- function(Y_t, J, T_max, early_null_periods = 0) {
  return(lapply(1:J,
                function(j) generate_Ps_j(Y_t, T_max, early_null_periods)))
}

#' @title Generate transition matrices for a transition pattern
#' @description Generate transition matrices for a transition pattern
#'
#' @param Y_t An outcome space object
#' @param T_max Maximum number of time periods in which transition matrices exist
#' @param early_null_periods Number of null periods before transition matrix appears; default = 0
#' matrix and initial outcome dists should be bounded below by some positive number.
#' Default value is FALSE; when TRUE, every element of transition matrices and
#' initial outcome dists are guaranteed to be strictly positive.
#'
#' @return A (T_max+early_null_periods-1)-length list of transition matrices for each time period
generate_Ps_j <- function(Y_t, T_max, early_null_periods = 0) {
  # Generate
  Ps_control <- lapply(rep(list(Y_t$names),(T_max-1)), generate_transition_matrix)

  # Add null periods if necessary
  if (early_null_periods > 0) {
    Ps_control <- c(vector("list", length = early_null_periods),
                    Ps_control)
  }

  # Return list
  return(Ps_control)
}

DECIMAL_POINT_MINIMAL <- 0.05
get_sample_size_for_model <- function() floor(1/DECIMAL_POINT_MINIMAL)
sample_for_model <- function(values) {
  sample(values, get_sample_size_for_model(), replace = TRUE)
}
round_to_decimal_point <- function(x) {
  return (round(x * 1/DECIMAL_POINT_MINIMAL) * DECIMAL_POINT_MINIMAL)
}

#' Generate a Random Transition Matrix from Specified Dimensions
#'
#' This function generates a random transition matrix of dimensions
#' \code{length(outcomes_past)} x \code{length(outcomes_future)},
#' where each row is a valid probability
#' mass function (PMF). Each row sums exactly to 1.
#'
#' @param outcomes_past A vector specifying the names of past states (rows).
#' @param outcomes_future A vector specifying the names of future states (columns).
#' If NULL, defaults to \code{K_past}.
#'
#' @return A numeric matrix of dimension \code{K_past} by \code{K_future},
#'   where each row represents a randomly generated PMF.
#'
#' @examples
#' transition_matrix <- generate_transition_matrix_from_K(3)
#' print(transition_matrix)
generate_transition_matrix <- function(outcomes_past, outcomes_future = NULL) {
  # If K_past is not provided, assume it is the same as K_future
  if (is.null(outcomes_future)) {
    outcomes_future <- outcomes_past
  }

  # Generate a random PMF for each element in the transition matrix
  P <- t(sapply(rep(list(outcomes_future), length(outcomes_past)), generate_pmf))
  rownames(P) <- outcomes_past
  P
}

#' Generate a Random Probability Mass Function (PMF) for an Arbitrary K
#'
#' This function generates a random probability mass function (PMF) vector of length \code{length(outcomes)}
#' where each element is a multiple of a specified granularity that ensures the sum of the vector is exactly 1.
#' The granularity is automatically determined based on the value of \code{length(outcomes)}.
#'
#' @param outcomes A vector specifying the names of states.
#' @return A numeric vector of length \code{K} representing a PMF where each element is a multiple of the chosen granularity,
#' and the sum of all elements equals 1.
#' @examples
#' pmf <- generate_pmf_from_K(4)
#' print(pmf)
generate_pmf <- function(outcomes) {
  K <- length(outcomes)
  granularity <- find_granularity(K)

  # If K = 1, return a vector of 1
  if (K == 1) {
    return(1)
  }

  # Initialize the vector with zeros
  pmf <- rep(0, K)

  # Randomly distribute units of the chosen granularity across the K elements
  remaining_mass <- 1

  for (i in 1:(K-1)) {
    max_mass <- min(remaining_mass - (K - i) * granularity, 1 - granularity)

    pmf[i] <- sample(seq(granularity, max(granularity, max_mass), by = granularity), 1)
    remaining_mass <- remaining_mass - pmf[i]
  }

  # Assign the remaining mass to the last element
  pmf[K] <- remaining_mass
  names(pmf) <- outcomes

  return(pmf)
}

#' Determine the Granularity for the PMF Based on the Length of the Vector
#'
#' This function calculates the granularity \code{10^{-J}5} such that \code{K * granularity <= 1}.
#' If K <= 2, use granularity of 0.1. Otherwise, calculate J such that \code{K * 10^{-J}5 <= 1}.
#' The granularity ensures that the generated PMF vector elements are multiples of this value,
#' and their sum is exactly 1.
#'
#' @param K An integer representing the number of categories in the PMF.
#' @return A numeric value representing the granularity to be used for generating the PMF.
#' @examples
#' granularity <- find_granularity(4)
#' print(granularity)
find_granularity <- function(K) {
  J <- ceiling((log(K) + log(5))/log(10))
  return(min(1 * 10^(-J), 0.1))
}
