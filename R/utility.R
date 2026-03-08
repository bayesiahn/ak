#' Get Option Value from Options List
#'
#' Generic helper function to retrieve an option value from an options list.
#' Returns the specified default if the option is not found or opts is not a list.
#'
#' @param opts An optional list containing algorithm options.
#' @param key The name of the option to retrieve.
#' @param default The default value to return if the option is not found.
#'
#' @return The option value if found, otherwise the default value.
#'
#' @examples
#' opts <- list(max_iter = 200, tol = 1e-5)
#' get_option(opts, "max_iter", 100)  # Returns 200
#' get_option(opts, "missing_key", 50)  # Returns 50
#' get_option(NULL, "max_iter", 100)  # Returns 100
#'
#' @export
get_option <- function(opts, key, default) {
  if (is.list(opts) && key %in% names(opts)) {
    return(opts[[key]])
  }
  return(default)
}

#' Check if a Value is Valid (Non-NULL and Non-NA)
#'
#' Validates if a value is neither NULL nor NA. This is useful for checking
#' whether a value has been properly set before using it.
#'
#' @param value The value to be validated.
#' @return A logical value: TRUE if the value is valid, FALSE otherwise.
#'
#' @examples
#' is_valid_value(5)       # Returns TRUE
#' is_valid_value(NULL)    # Returns FALSE
#' is_valid_value(NA)      # Returns FALSE
#' is_valid_value("")      # Returns TRUE (empty string is valid)
#'
#' @export
is_valid_value <- function(value) {
  !is.null(value) && !anyNA(value)
}

#' Check if a String Variable Name is Valid
#'
#' Validates if a variable name is a non-empty string. This is useful for
#' ensuring that variable names are correctly specified in specifications.
#'
#' @param varname A character string representing the variable name to be
#'   validated.
#' @return A logical value: TRUE if valid (non-NULL, non-NA, non-empty string),
#'   FALSE otherwise.
#'
#' @examples
#' is_valid_varname("treatment")  # Returns TRUE
#' is_valid_varname("")           # Returns FALSE
#' is_valid_varname(NULL)         # Returns FALSE
#' is_valid_varname(NA)           # Returns FALSE
#'
#' @export
is_valid_varname <- function(varname) {
  !is.null(varname) && !is.na(varname) && nchar(varname) > 0
}

#' Convert a vector of outcomes to a vector of indices
#' @param y A vector of outcomes.
#' @param y_values A vector of possible outcomes.
#' @return A vector of indices.
#' @export
#' @examples
#' y_to_y_indices(c("a", "b", "c"), c("a", "b", "c"))
#' y_to_y_indices(c("a", "b", "c"), c("a", "b"))
#' y_to_y_indices(c("a", "b", "c"), c("a", "b", "c", "d"))
y_to_y_indices <- function(y, y_values) {
  if (length(y_values) < 1) {
    return(NA)
  }
  indices <- sapply(y, function(y_i) which(y_values == y_i))
  if (!is.null(dim(y))) {
    indices <- matrix(indices, nrow = dim(y)[1])
  }
  return(indices)
}

#' Estimate Transition Indicators for Each Unit and Time Period
#'
#' This function computes the transition indicators for each unit across
#' different time periods. It returns a T-length list where each element is an
#' N-list of K by K matrices. Each matrix represents the transition pattern
#' of an individual unit in a particular period.
#'
#' @param y_indices_matrix A numeric matrix of N rows and T+1 columns,
#'        where each row represents a sequence of states (as indices) for
#'        an individual unit over time.
#' @param y_values A numeric vector representing the possible states.
#' @param current_period_to_be_fixed An integer representing the time period to be
#' fixed for current states; default is NULL
#'
#' @return A list of length T, where each element is a list of N transition matrices
#'         (each of size K by K). Each matrix in the list corresponds to a single
#'         unit's transition from one state to another in a given time period.
#'
#' @examples
#' y_indices_matrix <- matrix(c(1, 2, 2, 1, 2, 3), nrow = 2)
#' y_values <- 1:3
#' transitions <- transitions_to_transition_indicators(y_indices_matrix, y_values)
#' transitions
#'
#' @export
transitions_to_transition_indicators <- function(y_indices_matrix, y_values, current_period_to_be_fixed = NULL) {
  # Number of units and time periods
  N <- nrow(y_indices_matrix)
  T <- ncol(y_indices_matrix) - 1
  K <- length(y_values)

  # Return empty list if there are no units or more than two time periods
  if (N == 0 || T <= 0) {
    return(list())
  }

  # Initialize the list of lists for transition matrices
  transition_list <- vector("list", T)

  # Loop through each time period
  for (t in 1:T) {
    # Extract the current and next state indices for all units
    current_period <- ifelse(is.null(current_period_to_be_fixed), t, current_period_to_be_fixed)
    current_states <- y_indices_matrix[, current_period]
    next_states <- y_indices_matrix[, t + 1]

    # Vectorized operation to create transition matrices for all units
    Ps_control_t <- lapply(1:N, function(i) {
      # Create a K by K transition matrix with all zeros
      transition_matrix <- matrix(0, nrow = K, ncol = K)

      # Set the transition from current to next state for unit i
      transition_matrix[current_states[i], next_states[i]] <- 1

      return(transition_matrix)
    })

    # Assign the list of matrices for time t to the main list
    transition_list[[t]] <- Ps_control_t
  }

  return(transition_list)
}

#' Create Initial Outcome Indicators for Each Unit
#'
#' This function generates initial outcome indicators for each unit in the
#' `y_indices_matrix`, using a vectorized approach for efficiency. It returns
#' an N by K matrix, where each row represents the initial state of each unit
#' as a one-hot encoded vector.
#'
#' @param y_indices_matrix A numeric matrix of N rows and T+1 columns, where
#'        each row represents a sequence of states (as indices) for an individual
#'        unit over time. The first column should represent the initial state.
#' @param y_values A numeric vector representing the possible states.
#'
#' @return An N by K matrix, where each row is a one-hot encoded representation
#'         of the initial state of the corresponding unit.
#'
#' @examples
#' y_indices_matrix <- matrix(c(1, 2, 2, 3, 1, 2), nrow = 2)
#' y_values <- 1:3
#' initial_outcome_indicators <- y_indices_to_initial_outcome_indicators(y_indices_matrix, y_values)
#' initial_outcome_indicators
#'
#' @export
y_indices_to_initial_outcome_indicators <- function(y_indices_matrix, y_values) {
  N <- nrow(y_indices_matrix)
  K <- length(y_values)

  # Create a matrix of zeros with dimensions N x K
  initial_outcome_matrix <- matrix(0, nrow = N, ncol = K)

  # Set the appropriate indices to 1
  if (ncol(y_indices_matrix) > 0) {
    initial_states <- y_indices_matrix[, 1]
    initial_outcome_matrix[cbind(1:N, initial_states)] <- 1
  }

  return(initial_outcome_matrix)
}

#' Multiply Matrices in a List
#'
#' This function takes a list of matrices and returns the product of these matrices.
#' It assumes that the matrices are in a sequence that allows for matrix multiplication
#' (i.e., the number of columns in one matrix equals the number of rows in the next).
#'
#' @param matrix_list A list of matrices.
#'
#' @return The product of the matrices in the list.
#'
#' @examples
#' mat1 <- matrix(1:4, nrow = 2, ncol = 2)
#' mat2 <- matrix(1:4, nrow = 2, ncol = 2)
#' result <- multiply_matrices_in_list(list(mat1, mat2))
#'
#' @export
multiply_matrices_in_list <- function(matrix_list) {
  Reduce(`%*%`, matrix_list)
}


#' Concatenate Transition Matrices Lists
#'
#' This function concatenates two lists of transition matrices (Ps). It assumes that
#' both lists have the same length and correspond to the same number of transition patterns.
#'
#' @param Ps_pre A list of transition matrices before a certain time point or event.
#' @param Ps_post A list of transition matrices after a certain time point or event.
#'
#' @return A list of concatenated transition matrices for each transition pattern.
#'
#' @examples
#' Ps_pre <- list(matrix(1:4, nrow = 2), matrix(5:8, nrow = 2))
#' Ps_post <- list(matrix(9:12, nrow = 2), matrix(13:16, nrow = 2))
#' Ps <- concatenate_Ps(Ps_pre, Ps_post)
#'
#' @export
concatenate_Ps <- function(Ps_pre, Ps_post) {
  # Return the pre-transition matrices if post-transition matrices are NULL
  if (is.null(Ps_post)) {
    return(Ps_pre)
  }

  # Determine the number of transition patterns (J) from the pre-treatment transition matrices
  J <- length(Ps_pre)

  # Check if the number of transition patterns in pre- and post-treatment transition matrices are the same
  if (J != length(Ps_post)) {
    stop("Two lists of Ps must have the same length.")
  }

  # Initialize a list to store the concatenated transition matrices
  Ps <- vector("list", length = J)

  # Concatenate the transition matrices for each transition pattern
  for (j in 1:J) {
    # For each transition pattern, concatenate the corresponding pre- and post-treatment matrices
    Ps[[j]] <- c(Ps_pre[[j]], Ps_post[[j]])
  }

  # Return the concatenated list of transition matrices
  return(Ps)
}


#' Extract y matrix and g vector from long panel data
#' @param data long panel data
#' @param yname name of dependent variable
#' @param gname name of group membership variable
#' @param idname name of id variable
#' @param tname name of time variable
#' @param cluster_varname name of cluster variable (optional); default is NULL
#' @return list with y matrix and g vector
#' @examples
#' data <- data.frame(id = c(1,1,1,2,2,2), t = c(1,2,3,1,2,3),
#' y = c(1,2,3,4,5,6), g = c(1,1,1,2,2,2))
#' long_df_to_y_and_g(data, "y", "g", "id", "t")
#' @export
long_df_to_y_and_g <- function(data, yname, gname, idname, tname, cluster_varname = NULL) {
  # Extract wide dfs
  wide_df_y <- long_df_to_wide(data, yname, idname, tname)
  wide_df_g <- long_df_to_wide(data, gname, idname, tname)

  # Extract y
  y <- as.matrix(wide_df_y)

  # Extract g
  g_unscaled <- apply(wide_df_g, 1, dplyr::first)
  t_range_unscaled <- sort(unique(dplyr::pull(data, tname)))
  g_range_unscaled <- sort(unique(c(g_unscaled)))

  # Extract c if needed
  c <- NULL
  if (!is.null(cluster_varname)) {
    wide_df_c <- long_df_to_wide(data, cluster_varname, idname, tname)
    # Choose only first column as clusters are invariant over time
    c <- wide_df_c[, 1]
  }

  # If g does not appear in the time range, assign zeros
  g_range <- match(g_range_unscaled, t_range_unscaled)
  g_range[is.na(g_range)] <- 0
  g <- g_range[match(c(g_unscaled), g_range_unscaled)]

  return(list(y = y, g = g, c = c))
}

#' Change long panel data to wide
#' @param data long panel data
#' @param varname name of variable to extract
#' @param idname name of id variable
#' @param tname name of time variable
#' @return Matrix of variable varname
#' @examples
#' data <- data.frame(id = c(1,1,1,2,2,2), t = c(1,2,3,1,2,3),
#' y = c(1,2,3,4,5,6), g = c(1,1,1,2,2,2))
#' data_to_y_and_gs(data, "y", "g", "id", "t")
long_df_to_wide <- function(data, varname, idname, tname) {
  wide_df <- as.data.frame(data) %>%
    dplyr::select(dplyr::all_of(varname),
                  dplyr::all_of(idname),
                  dplyr::all_of(tname)) %>%
    stats::reshape(idvar = idname, timevar = tname, direction = "wide")
  outcome_cols <- grep(varname, colnames(wide_df), value = TRUE)
  outcome_matrix <- as.matrix(wide_df[, outcome_cols])
  outcome_matrix
}
