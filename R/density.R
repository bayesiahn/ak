#' Compute Probability Mass Function
#'
#' This function computes the probability mass function (PMF) of the given outcomes `y_t` based on the outcome space `Y_t`.
#'
#' @param y_t Numeric vector representing outcomes for each unit.
#' @param Y_t A list containing 'values', a sorted numeric vector, and 'type', a string indicating the outcome type ('discrete' or 'continuous').
#' @return A numeric vector representing the probability mass function (PMF).
#' @examples
#' y_t <- c(1, 1, 2, 2, 2, 3)
#' Y_t <- list(values = c(1, 2, 3), type = 'discrete')
#' result <- get_pmf(y_t, Y_t)
#' # Should return c(1/6, 3/6, 2/6) for values 1, 2, and 3 respectively
#' @export
get_pmf <- function(y_t, Y_t) {
  N <- length(y_t)
  Y_t_indices <- 1:length(Y_t$values)
  y_t_by_Y_t_index <- get_y_t_by_Y_t_index(y_t, Y_t)
  pmf <- sapply(Y_t_indices, function(i) sum(y_t_by_Y_t_index == i)) / N
  names(pmf) <- get_Y_t_names(Y_t)
  return(pmf)
}
