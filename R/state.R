#' Get the Default Grid Length
#'
#' This function computes the default grid length given all the observations
#' in the outcome space \( Y_t \). The grid length is calculated as the ceiling
#' of the ratio of the total number of observations to a specified ratio, which defaults to 20.
#'
#' @param Y_t_all A numeric vector containing all the observations in the outcome space \( Y_t \).
#' @param ratio An optional numeric value to divide the total number of observations. Defaults to 20.
#'
#' @return An integer representing the default grid length.
#'
#' @examples
#' Y_t_all <- c(1, 2, 3, 4, 5, 6)
#' get_default_grid_length(Y_t_all) # Should return 1
#' get_default_grid_length(Y_t_all, ratio = 6) # Should return 1
#'
#' @export
get_default_grid_length <- function(Y_t_all, ratio = 10) {
  ceiling(length(Y_t_all) / ratio)
}

#' Check if Outcome Space \( y_t \) is Discrete
#'
#' This function checks whether the given outcome space \( y_t \) is discrete.
#' Specifically, it returns TRUE if there are duplicate values in \( y_t \),
#' otherwise FALSE.
#'
#' @param y_t A numeric vector representing the outcome space \( y_t \).
#'
#' @return A logical value indicating whether \( y_t \) is discrete or not.
#'
#' @examples
#' is_discrete(c(1, 1, 2)) # Should return TRUE
#' is_discrete(c(1, 2, 3)) # Should return FALSE
#' @export
is_discrete <- function(y_t) {
  y_t_unique <- unique(y_t)
  (length(y_t_unique) != length(y_t)) && length(y_t_unique) < 10
}
#' Generate Outcome Space Details
#'
#' This function takes a vector of outcomes `y_t` and optionally, a list of parameters,
#' to generate a list describing the outcome space details, which include the type ('discrete' or 'continuous')
#' and a vector of unique values or quantiles.
#'
#' @param outcomes_observed Numeric / string vector of outcomes for each unit.
#' @param outcomes_of_interest String vector of outcomes of interest; default is the last element of sort(unique(y_t)) 
#' @return A list containing 'values', a sorted numeric vector of unique values or quantiles,
#'         and 'type', a character string indicating whether the outcome space is 'discrete' or 'continuous'.
#' @examples
#' y_t_discrete <- c(1, 2, 2, 3)
#' result_discrete <- generate_Y_t(y_t_discrete)
#' print(result_discrete$type)  # Should print 'discrete'
#' print(result_discrete$values)  # Should print c(0, 0, 1)
#' result_discrete <- generate_Y_t(y_t_discrete, 2)
#' print(result_discrete$values)  # Should print c(0, 1, 0)
#' @export
generate_Y_t <- function(outcomes_observed, outcomes_of_interest = NULL) {
  names <- sort(unique(outcomes_observed))

  type <- "discrete"
  lower_bound <- min(names)

  # Assign values; choose the last member as one and rest as zero
  if (is.null(outcomes_of_interest)) {
    outcomes_of_interest <- names[length(names)]
  }
  values <- ifelse(names %in% outcomes_of_interest, 1, 0)

  list(values = values, names = names, type = type, lower_bound = lower_bound)
}

#' Generate names for a set of intervals
#'
#' Given a vector of interval endpoints, this function returns a character vector of interval names in the form "[a, b)".
#'
#' @param interval_vector A numeric vector of sorted interval endpoints.
#' @return A character vector of interval names.
#'
#' @examples
#' get_interval_names(c(0, 0.1, 0.3, 0.5))
#' # should return c("[0.00, 0.10)", "[0.10, 0.30)", "[0.30, 0.50]")
#' @export
get_interval_names <- function(interval_vector) {
  n_intervals <- length(interval_vector) - 1
  if (n_intervals <= 0) {
    return(character(0))
  }

  interval_names <- character(n_intervals)

  for (i in 1:n_intervals) {
    lower_bound <- sprintf("%.2f", interval_vector[i])
    upper_bound <- sprintf("%.2f", interval_vector[i + 1])
    interval_names[i] <- paste0("[", lower_bound, ", ",
                                upper_bound, ")")
    if (i == n_intervals) {
      interval_names[i] <- paste0("[", lower_bound, ", ",
                                  upper_bound, "]")
    }
    if (lower_bound == upper_bound) {
      interval_names[i] <- lower_bound
    }
  }

  return(interval_names)
}

#' Get Lower Bound of Outcome Space
#'
#' Extracts the lower bound from an outcome space object.
#'
#' @param Y_t An outcome space object with optional `lower_bound` field.
#' @return The lower bound value.
#' @keywords internal
get_lower_bound <- function(Y_t) {
  if ("lower_bound" %in% names(Y_t)) {
    return(Y_t$lower_bound)
  } else {
    return(min(Y_t$values))
  }
}

#' Check if Outcome Space is Discrete
#'
#' @param Y_t An outcome space object with a `type` field.
#' @return TRUE if discrete, FALSE otherwise.
#' @keywords internal
is_Y_t_discrete <- function(Y_t) {
  return(Y_t$type != "continuous")
}

#' Check if Outcome Space is Continuous
#'
#' @param Y_t An outcome space object with a `type` field.
#' @return TRUE if continuous, FALSE otherwise.
#' @keywords internal
is_Y_t_continuous <- function(Y_t) {
  return(Y_t$type == "continuous")
}

#' Generate names for Y_t
#'
#' Given a list Y_t that describes the type and values of an outcome,
#' this function generates names for the outcomes.
#' The function differentiates between continuous and categorical types.
#'
#' @param Y_t A list containing the type ("continuous" or otherwise) and values for Y_t.
#' For continuous type, it should also contain the lower bound.
#' @return A character vector or numeric vector, depending on the type in Y_t.
#'
#' @examples
#' get_Y_t_names(list(type = "continuous", lower_bound = 0, values = c(0.1, 0.2)))
#' # should return c("[0.00, 0.10)", "[0.10, 0.20)")
#'
#' get_Y_t_names(list(type = "categorical", values = c("A", "B")))
#' # should return c("A", "B")
#' @export
get_Y_t_names <- function(Y_t) {
  Y_t$names
}

#' Get Conditional State Matrix
#'
#' This function returns a matrix where each element is a concatenation of the current and previous `lags` elements
#' from the input matrix `y`. If `lags = 1`, each element will be a concatenation of the current and previous element.
#' If `lags = 2`, each element will be a concatenation of the current and previous two elements, and so on.
#'
#' @param y A numeric matrix.
#' @param lags An integer specifying the number of lags. Default is 1.
#'
#' @return A matrix of the same dimensions as `y` where each element is a string of concatenated values
#' separated by underscores.
#'
#' @examples
#' y <- matrix(1:6, nrow = 2, byrow = TRUE)
#' get_conditional_state_matrix(y, lags = 1)
#' # Returns:
#' #      [,1]   [,2]   [,3]
#' # [1,] NA    "1_2" "2_3"
#' # [2,] NA    "4_5" "5_6"
#'
#' get_conditional_state_matrix(y, lags = 2)
#' # Returns:
#' #      [,1]   [,2]   [,3]
#' # [1,] NA    NA    "1_2_3"
#' # [2,] NA    NA    "4_5_6"
#'
#' @export
get_conditional_state_matrix <- function(y, lags = 1) {
  if (lags == 0) {
    return(y)
  }
  nrow_y_indices_matrix <- nrow(y)
  ncol_y_indices_matrix <- ncol(y)
  y_conditional_state_matrix <- matrix(NA,
                                       nrow = nrow_y_indices_matrix, ncol = ncol_y_indices_matrix)
  # Vectorize across rows for each time point (avoids nested loop)
  # Only process if there are enough columns for the given lags
  if (lags + 1 <= ncol_y_indices_matrix) {
    for (t in (lags + 1):ncol_y_indices_matrix) {
      y_conditional_state_matrix[, t] <- apply(y[, (t - lags):t, drop = FALSE], 1,
                                               paste, collapse = "_")
    }
  }
  return(y_conditional_state_matrix)
}

