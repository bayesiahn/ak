#' Get Number of Latent Groups
#'
#' Extracts the number of latent groups (J) from a transition model.
#'
#' @param transition_model A transition model object containing priors.
#' @return Integer number of latent groups.
#' @keywords internal
get_J <- function(transition_model) length(transition_model$priors)

#' Get Maximum Time Period
#'
#' Extracts the maximum number of time periods (T) from a transition model.
#'
#' @param transition_model A transition model object containing transition matrices.
#' @return Integer number of time periods.
#' @keywords internal
get_T_max <- function(transition_model) length(transition_model$Ps_control[[1]])

#' Generate a sample from a transition model
#'
#' @param transition_model A list representing the transition model. It should include the following elements:
#'   - Y: An outcome space object.
#'   - priors: A numeric vector representing the prior probabilities for each group.
#'   - Ps_control: A list of transition matrices, one for each group.
#'   - initial_outcome_dist: A list of initial outcome distributions, one for each group.
#' @param N The number of individuals to sample. Default is 20.
#' @return A list with the following elements:
#'   - y: A matrix of outcomes for all individuals.
#'   - y_its: A list of matrices of outcomes, one for each group.
#'   - N_js: A table showing the number of individuals in each group.
#'   - transition_model: The input transition model.
#' @export
generate_sample <- function(transition_model, N = 20) {
  J <- get_J(transition_model)
  prob_treated <- transition_model$prob_treated

  # If treated_period <= 0, then no one is treated
  if (transition_model$treated_period <= 0) {
    prob_treated <- 0
  }
  probs_treated_by_type <- compute_probs_treated_by_type(transition_model$priors_treated,
                                                   transition_model$priors,
                                                   prob_treated)


  # Draw and count number of individuals in each group
  N_js <- table(sample(factor(1:J), N, replace = TRUE, prob = transition_model$priors))

  # Draw and count number of individuals in each group who are treated
  # First row represents the number of individuals in each group who are not treated
  # Second row represents the number of individuals in each group who are treated
  N_gjs_table <- sapply(1:J, function(j)
    table(sample(factor(c('Control',paste0('Treated (period ', transition_model$treated_period, ')'))), N_js[j],
                 replace = TRUE, prob = c(1-probs_treated_by_type[j], probs_treated_by_type[j]))))

  # Generate outcomes for each treatment group
  y_its_control <- lapply(1:J, function(j) {
    generate_sample_j(transition_model$Y, N_gjs_table[1,j],
                      transition_model$Ps_control[[j]],
                      transition_model$pmfs_initial_control[[j]])})
  y_its_treated <- lapply(1:J, function(j) {
    generate_sample_j(transition_model$Y, N_gjs_table[2,j],
                      transition_model$Ps_treated[[j]],
                      transition_model$pmfs_initial_treated[[j]])})

  # Collect them as a matrix
  y_control <- do.call(rbind, y_its_control)
  y_treated <- do.call(rbind, y_its_treated)
  y <- y_control
  g <- rep(0, N)
  j <- unlist(sapply(1:J, function(j_index) rep(j_index, N_gjs_table[1,j_index])))
  if (length(y_treated) > 0) {
    y <- rbind(y_control, y_treated)
    g <- c(rep(0, nrow(y_control)), rep(transition_model$treated_period, nrow(y_treated)))
    j <- c(j, unlist(sapply(1:J, function(j_index) rep(j_index, N_gjs_table[2,j_index]))))
  }
  y <- as.matrix(y, nrow = N)
  # Return outcomes
  list(y = y, g = g, j = j,
       y_its_control = y_its_control, y_its_treated = y_its_treated,
       N_js = N_js, N_gjs_table = N_gjs_table,
       transition_model = transition_model)
}

#' Generate a sample for a specific group
#'
#' @param Y An outcome space object.
#' @param N_j The number of individuals to sample for this group.
#' @param Ps_j T-list of transition matrices for this group.
#' @param initial_outcome_dist_j J-length vector of initial outcome distribution for this group.
#' @return A matrix of outcomes for the individuals in this group.
#' @export
generate_sample_j <- function(Y, N_j,
                              Ps_j, initial_outcome_dist_j) {
  T_max <- length(Ps_j) + 1
  if (N_j <= 0) {
    return(matrix(numeric(0), ncol = T_max))
  }
  # Generate initial outcomes
  y_0 <- sample(Y$names, N_j, replace = TRUE, prob = initial_outcome_dist_j)
  y_its <- matrix(y_0, ncol = 1)

  # Generate outcomes
  for (t in 1:(T_max-1)) {
    # Get indices of past outcomes
    y_past_indices <- y_to_y_indices(y_its[, t], Y$names)

    # Sample outcomes
    P_jt <- Ps_j[[t]]
    y_t <- sapply(y_past_indices, function(y_past_index) {
                  sample(Y$names, 1, replace = TRUE,
                  prob = P_jt[y_past_index, ])})

    # Add outcomes to matrix
    y_its <- cbind(y_its, y_t)
  }

  # Return outcomes; if no character, convert to numbers
  y_its
}
