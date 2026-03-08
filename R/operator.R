#' Calculate the midpoints between intervals
#'
#' @param intervals Numeric vector containing sorted intervals.
#' @return Returns a numeric vector containing the midpoints between intervals.
#'
#' @examples
#' get_midpoints(c(0.1, 0.3, 0.5))  # should return c(0.2, 0.4)
#' @export
get_midpoints <- function(intervals) {
  if (!is.numeric(intervals) || length(intervals) < 2) {
    stop("Invalid input: provide a numeric vector of length at least 2.")
  }
  mid_values <- (intervals[-1] + intervals[-length(intervals)]) / 2
  return(mid_values)
}


#' Calculate the expectation operator \( h_t \)
#'
#' This function computes the expectation operator \( h_t \), a column vector
#' such that t(f) %*% h_t gives the expected value.
#'
#' @param Y_t A list containing the 'values' vector and 'type' string.
#'           'values' is a sorted numeric vector and 'type' indicates whether the
#'           outcome space is 'discrete' or 'continuous'.
#'
#' @return A numeric vector \( h_t \) serving as the expectation operator.
#'
#' @examples
#' Y_t_discrete <- list(values = c(1, 2, 3), type = 'discrete')
#' get_expectation_operator(Y_t_discrete)
#'
#' Y_t_continuous <- list(values = c(0.1, 0.3, 0.5), type = 'continuous')
#' get_expectation_operator(Y_t_continuous)
#' @export
get_expectation_operator <- function(Y_t) {
  # Ensure 'Y_t' is a list and contains 'values' and 'type'
  if (!is.list(Y_t) || !all(c("values", "type") %in% names(Y_t))) {
    stop("Invalid Y_t")
  }

  values <- Y_t$values

  # Depending on 'type', compute h_t differently
  if (Y_t$type == "discrete") {
    h_t <- values
  } else if (Y_t$type == "continuous") {
    h_t <- get_midpoints(c(Y_t$lower_bound, values))
  } else {
    stop("Unknown type specified for Y_t")
  }

  return(h_t)
}
