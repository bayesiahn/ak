

#' Convert Priors for Treated to Probabilities of Being Treated for Each Latent Group
#'
#' This function computes the probabilities of being treated based on the priors for treated individuals,
#' overall priors, and the overall probability of being treated. It ensures that the computed probabilities
#' do not exceed 1.
#'
#' @param priors_treated A vector of prior probabilities for each group given the unit is treated.
#' @param priors A vector of overall prior probabilities for each group.
#' @param prob_treated The overall probability of being treated in the sample.
#'
#' @return A vector of probabilities of being treated for each group.
#'
#' @examples
#' priors_treated <- c(0.6, 0.4)
#' priors <- c(0.5, 0.5)
#' prob_treated <- 0.8
#' probs_treated_by_type <- compute_probs_treated_by_type(priors_treated, priors, prob_treated)
compute_probs_treated_by_type <- function(priors_treated, priors, prob_treated) {
  if (prob_treated <= 0) {
    priors_treated <- rep(0, length(priors))
  }
  probs_treated_by_type <- priors_treated * prob_treated / priors
  if (any(probs_treated_by_type > 1)) {
    stop("The probability of being treated is greater than 1 for at least one group; \n
         check priors_treated (P(Z_i = j | D_i = 1)) and priors (P(Z_i = j)) so that \n
         priors_treated[j] * prob_treated / priors[j] <= 1 for all j. ")
  }
  return(probs_treated_by_type)
}


#' Get Fields Indexed by Latent Group (J)
#'
#' Returns the list of TransitionModel field names that are indexed by latent
#' group J and need reordering when model groups are reordered.
#'
#' @return A character vector of field names.
#' @keywords internal
get_j_indexed_fields <- function() {
  c(
    # Core priors (always present)
    "priors",
    "priors_treated",
    # Transition matrices (always present)
    "Ps_control",
    "Ps_treated",
    "Ps_control_empirical",
    "Ps_treated_empirical",
    # Initial outcome distributions
    "pmfs_initial_control",
    "pmfs_initial_treated",
    # Pre-treatment outcome distributions
    "pmfs_pre_treatment_control",
    "pmfs_pre_treatment_treated",
    # Optional transition matrices from different reference periods
    "Ps_control_from_pre",
    "Ps_treated_from_pre",
    "Ps_control_from_0",
    "Ps_treated_from_0"
  )
}

#' Get Matrix Fields with J Columns
#'
#' Returns the list of TransitionModel field names that are matrices with
#' columns indexed by latent group J (e.g., posteriors, weights).
#'
#' @return A character vector of field names.
#' @keywords internal
get_j_column_fields <- function() {
  c("posteriors", "weights")
}

#' Reorder Transition Model Estimate by Priors
#'
#' Reorders the components of a TransitionModel object based on the priors,
#' aligning all elements of the transition model according to the order of
#' the priors.
#'
#' @param transition_model A TransitionModel object to be reordered.
#' @return A reordered TransitionModel object.
#' @export
reorder_transition_model <- function(transition_model) {
  if (!inherits(transition_model, "TransitionModel")) {
    stop("Input must be a TransitionModel object.")
  }

  # Extract priors and determine the order
  priors <- transition_model$priors
  order_indices <- order(priors, decreasing = TRUE)

  # Reorder J-indexed fields (lists/vectors indexed by latent group)
  j_indexed_fields <- get_j_indexed_fields()
  for (field in j_indexed_fields) {
    if (!is.null(transition_model[[field]])) {
      transition_model[[field]] <- transition_model[[field]][order_indices]
    }
  }

  # Reorder matrix fields with J columns
  j_column_fields <- get_j_column_fields()
  for (field in j_column_fields) {
    if (!is.null(transition_model[[field]])) {
      transition_model[[field]] <- transition_model[[field]][, order_indices]
    }
  }

  return(transition_model)
}



#' Calculate Occurrence Rate for Control and Treated Groups
#'
#' Computes the occurrence rates for control and treated groups across all
#' groups defined by the priors. This function uses the treated and overall
#' priors along with the overall probability of being treated to calculate
#' the occurrence rate for both control and treated groups for each group `j`.
#'
#' @param p_y1d A J-length list of joint pmfs for (Y_1, D) represented by K by 2 matrices (1st column is for D = 0, 2nd column is for D = 1)
#' @param priors A vector of overall prior probabilities for each group.
#' @return A 2xJ matrix with occurrence rates for control (row 1) and treated (row 2) groups across all J groups.
#' @examples
#' p_y1d <- list(matrix(c(0.3, 0.2, 0.4, 0.1), ncol = 2),
#'               matrix(c(0.2, 0.1, 0.1, 0.6), ncol = 2))
#' priors <- c(0.5, 0.5)
#' occurance_rate <- get_occurence_rates(p_y1d, priors)
#' print(occurance_rate)
#' @export
get_occurence_rates <- function(p_y1d, priors) {
  # Retrieve the number of components
  J <- length(priors)

  # Retrieve probs_treated_by_type
  probs_not_treated_by_j <- sapply(p_y1d, function(p_y1d_j) sum(p_y1d_j[,1]))
  probs_treated_by_j <- sapply(p_y1d, function(p_y1d_j) sum(p_y1d_j[,2]))

  # Get the occurance rate
  occurance_rate <- matrix(0, nrow = 2, ncol = J)
  for (j in 1:J) {
    occurance_rate[1,j] <- probs_not_treated_by_j[j] * priors[j]
    occurance_rate[2,j] <- probs_treated_by_j[j] * priors[j]
  }

  # Attach colnames and rownames
  colnames(occurance_rate) <- 1:J
  rownames(occurance_rate) <- c("Control", "Treated")

  occurance_rate
}

#' Get Expected Outcomes Data Frame from Transition Model
#'
#' Computes expected outcomes for control and treated groups across all time periods
#' based on the specified transition model. It optionally aggregates these outcomes over
#' all groups. The function calculates both the immediate and future expected outcomes
#' using initial distributions and transition probabilities provided in the transition model.
#'
#' @param transition_model A list containing elements of the transition model, such as priors,
#'        transition matrices (`Ps_treated` and `Ps_control`), and outcome values (`Y$values`).
#'        The transition model should also include `treated_period`, the period in which treatment begins,
#'        and `T_max`, the maximum time period considered in the model.
#' @param aggregate Logical; if `TRUE`, the function aggregates outcomes over all groups (j),
#'        providing a summary outcome for each time period (`t`) and treatment status (`g`).
#'        If `FALSE`, outcomes are returned for each group and time period separately.
#'
#' @return A data frame containing columns `j` (group index), `t` (time period), `y` (expected outcome),
#'         `g` (treatment status), and `occurance_rate` (the rate of occurrence for each group and treatment status).
#'         The data frame presents outcomes for both control (g=0) and treated (g=`treated_period`) groups across
#'         all considered time periods.
#'
#' @examples
#' # Assuming a predefined transition_model
#' outcomes_df <- get_expected_outcomes_df(transition_model)
#' # To aggregate outcomes over all groups
#' aggregated_outcomes_df <- get_expected_outcomes_df(transition_model, aggregate = TRUE)
#'
#' @export
#' @importFrom dplyr left_join group_by reframe mutate
get_expected_outcomes_df <- function(transition_model, aggregate = FALSE) {
  # Get the weighted averages for treated units
  J <- length(transition_model$priors)
  g_max <- transition_model$treated_period
  T_max <- length(transition_model$Ps_control[[1]]) + 1
  Y_values <- transition_model$Y$values
  outcomes_control <- matrix(0, nrow = J, ncol = (T_max))
  outcomes_treated <- matrix(0, nrow = J, ncol = (T_max))

  # Occurence rates
  occurence_rates <- get_occurence_rates(transition_model$p_y1d, transition_model$priors)
  occurence_rates <- data.frame(occurance_rate = c(occurence_rates)) %>%
    cbind(expand.grid(g = c(0, g_max), j = 1:J))

  # Compute the difference in expected outcomes for each time period
  # Use cumulative products to avoid recomputation: O(T) instead of O(T²)
  for (j in 1:J) {
    outcomes_control[j,1] <- get_expected_outcome(transition_model$pmfs_initial_control[[j]],
                                                 Y_values)
    outcomes_treated[j,1] <- get_expected_outcome(transition_model$pmfs_initial_treated[[j]],
                                                  Y_values)

    # Initialize cumulative products for this latent group
    P_cumulative_treated <- NULL
    P_cumulative_control <- NULL

    for (t in 1:(T_max-1)) {
      # Update cumulative products incrementally
      if (t == 1) {
        P_cumulative_treated <- transition_model$Ps_treated_empirical[[j]][[1]]
        P_cumulative_control <- transition_model$Ps_control_empirical[[j]][[1]]
      } else {
        P_cumulative_treated <- P_cumulative_treated %*%
          transition_model$Ps_treated_empirical[[j]][[t]]
        P_cumulative_control <- P_cumulative_control %*%
          transition_model$Ps_control_empirical[[j]][[t]]
      }

      outcomes_control[j,t+1] <- get_expected_future_outcome(
        transition_model$pmfs_initial_control[[j]], P_cumulative_control, Y_values)
      outcomes_treated[j,t+1] <- get_expected_future_outcome(
        transition_model$pmfs_initial_treated[[j]], P_cumulative_treated, Y_values)
    }
  }

  # Create a data frame for the transition model
  jt_combinations <- expand.grid(j = 1:J, t = 1:T_max)
  outcomes_control_df <- jt_combinations %>%
    cbind(data.frame(y = as.vector(outcomes_control), g = 0))
  outcomes_treated_df <- jt_combinations %>%
    cbind(data.frame(y = as.vector(outcomes_treated), g = g_max))

  # Combine them together and merge with occurence rates
  outcomes_df <- rbind(outcomes_treated_df, outcomes_control_df) %>%
    dplyr::left_join(occurence_rates, by = c("j", "g"))

  # Get the weighted averages
  if (aggregate) {
    outcomes_df <- outcomes_df %>%
      dplyr::group_by(t, g) %>%
      dplyr::reframe(y = sum(occurance_rate*y)/sum(occurance_rate)) %>%
      dplyr::mutate(j = "Aggregate (All observations)")
  }

  outcomes_df
}

#' Summarize Transition Probabilities from Transition Model
#'
#' Generates a data frame that summarizes the transition probabilities for control
#' and treated groups for every possible combination of future and past outcomes
#' specified in the transition model. It provides insights into how likely transitions
#' between outcomes are under different conditions.
#'
#' @param transition_model A list containing the transition model elements, including
#'        transition matrices (`Ps_treated` and `Ps_control`), and a vector of
#'        outcome names (`Y$names`). The transition model should ideally include
#'        additional elements like `treated_period`, `T_max`, and `priors` to
#'        facilitate comprehensive analysis.
#' @param aggregate Logical; if `TRUE`, the function aggregates transition
#'        probabilities over all groups, providing a summary for each outcome
#'        transition. If `FALSE`, transition probabilities are returned for each
#'        group separately.
#'
#' @return A data frame containing columns for the past outcome (`past_outcome`),
#'         future outcome (`future_outcome`), transition probability (`probability`),
#'         treatment status (`treatment_status`), and optionally the group (`j`)
#'         if not aggregated. This summarizes the transition probabilities for all
#'         outcome combinations under both treated and control conditions.
#'
#' @examples
#' # Assuming a predefined transition_model with Y$names defined
#' transition_probs_df <- get_transition_probability_df(transition_model)
#' # To aggregate transition probabilities over all groups
#' aggregated_transition_probs_df <- get_transition_probability_df(transition_model, aggregate = TRUE)
#'
#' @export
#' @importFrom dplyr bind_rows group_by summarise mutate
get_transition_probability_df <- function(transition_model, aggregate = FALSE) {
  # Extract relevant elements from transition model
  Y_names <- transition_model$Y$names
  J <- length(transition_model$priors)
  g_max <- transition_model$treated_period

  # Occurence rates
  occurence_rates <- get_occurence_rates(transition_model$p_y1d, transition_model$priors)
  occurence_rates <- data.frame(occurance_rate = c(occurence_rates)) %>%
    cbind(expand.grid(g = c(0, g_max), j = 1:J))

  # Initialize empty list to store dfs
  dfs_list <- list()

  # Loop through control and treated to compute probabilities
  for (g in c(0, g_max)) {
    P_matrices <- transition_model$Ps_control_empirical
    if (g > 0) {
      P_matrices <- transition_model$Ps_treated_empirical
    }

    for (j in 1:J) {
      P_j <- P_matrices[[j]]
      for (t in 1:length(P_j)) {
        P_t <- P_j[[t]]
        df <- expand.grid(y_past = Y_names, y_current = Y_names)
        df$probability <- as.vector(P_t)
        df$j <- j
        df$g <- g
        df$t <- t + 1
        dfs_list[[length(dfs_list) + 1]] <- df
      }
    }
  }

  # Combine all data frames
  combined_df <- dplyr::bind_rows(dfs_list) %>%
    dplyr::left_join(occurence_rates, by = c("j", "g"))

  if (aggregate) {
    combined_df <- combined_df %>%
      dplyr::group_by(y_past, y_current, t, g) %>%
      dplyr::reframe(probability = sum(occurance_rate*probability)/sum(occurance_rate)) %>%
      dplyr::mutate(j = "Aggregate (All observations)")
  }

  combined_df
}

#' Modify Outcomes of Interest in a Model
#'
#' This function updates the outcome variables in a model object by
#' generating new outcome values based on the specified outcomes of interest.
#'
#' @param model A model object that contains outcome variable information.
#' @param outcomes_of_interest A vector of outcome names to update in the model.
#'
#' @return A modified model object with updated outcomes.
#'
#' @examples
#' model <- list(Y = list(names = c("outcome1", "outcome2"), values = c(0, 1)))
#' new_model <- modify_outcomes_of_interest(model, c("outcome1"))
#'
#' @export
modify_outcomes_of_interest <- function(model, outcomes_of_interest) {
  model$Y <- generate_Y_t(model$Y$names, outcomes_of_interest)

  # Update Y space in bootstrap estimates if available
  if (!is.null(model$bootstrap_estimates)) {
    model$bootstrap_estimates <- lapply(model$bootstrap_estimates,
                                        modify_outcomes_of_interest, outcomes_of_interest)
  }

  model
}

