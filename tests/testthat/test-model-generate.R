
test_that("generate_transition_model returns correct output", {
  # Define inputs
  outcomes <- c(0,1)
  J <- 2
  T_max <- 10

  # Call the function
  result <- generate_transition_model(outcomes, J, T_max)

  # Check that the result is a list
  expect_type(result, "list")

  # Check that the list has the correct names
  expect_equal(names(result), c("Y", "priors", "prob_treated",
                                "priors_treated", "treated_period",
                                "Ps_control", "Ps_treated",
                                "Ps_control_empirical", "Ps_treated_empirical",
                                "pmfs_initial_control", "pmfs_initial_treated",
                                "pmfs_pre_treatment_control",
                                "pmfs_pre_treatment_treated", "p_y1d"))

  # Check that Y is an outcome space object
  expect_type(result$Y, "list")

  # Check that priors is a numeric vector with J elements
  expect_type(result$priors, "double")
  expect_equal(length(result$priors), J)

  # Check that priors_treated is a numeric vector with J elements
  expect_type(result$priors_treated, "double")
  expect_equal(length(result$priors_treated), J)

  # Check that Ps_control is a list with J elements
  expect_type(result$Ps_control, "list")
  expect_equal(length(result$Ps_control), J)

  # Check that pmfs_initial_control is a list with J elements
  expect_type(result$pmfs_initial_control, "list")
  expect_equal(length(result$pmfs_initial_control), J)
})

test_that("pmfs_pre_treatment_control and _treated return correct output", {
  # Define inputs
  outcomes <- c(0,1)
  J <- 2
  T_max <- 10
  treated_period <- 3
  result <- generate_transition_model(outcomes, J, T_max, treated_period = treated_period)

  pmfs_pre_treatment_control <- result$pmfs_pre_treatment_control
  pmfs_pre_treatment_treated <- result$pmfs_pre_treatment_treated

  # Check that pmfs_pre_treatment_control _treated is a list with J elements
  expect_type(pmfs_pre_treatment_control, "list")
  expect_type(pmfs_pre_treatment_treated, "list")
  expect_equal(length(pmfs_pre_treatment_control), J)
  expect_equal(length(pmfs_pre_treatment_treated), J)

  # Check that the pmfs are equal
  for (j in 1:J) {
    pmf_pre_treatment_control_expected <- t(result$pmfs_initial_control[[j]]) %*%
      result$Ps_control[[j]][[1]] %>%
      as.numeric()
    pmf_pre_treatment_treated_expected <- t(result$pmfs_initial_treated[[j]]) %*%
      result$Ps_treated[[j]][[1]] %>%
      as.numeric()

    # Check that the pmfs are equal
    expect_equal(pmfs_pre_treatment_control[[j]], pmf_pre_treatment_control_expected)
    expect_equal(pmfs_pre_treatment_treated[[j]], pmf_pre_treatment_treated_expected)

    # Check if each jth element of p_y1d sums up to one
    expect_equal(sum(result$p_y1d[[j]]), 1)
  }
})

test_that("Check if p_y1d is valid", {
  Js <- 1:4
  outcomes_list <- list(c(0, 1), c(0, 1, 2), c(0, 1, 2, "A"))

  for (J in Js) {
    for (outcomes in outcomes_list) {
      result <- generate_transition_model(outcomes, J)

      for (j in 1:J) {
        # Check if each jth element of p_y1d sums up to one
        expect_equal(sum(result$p_y1d[[j]]), 1)
      }
    }
  }
})

test_that("generate_Ps returns correct output", {
  # Define inputs
  Y <- generate_Y_t(c(0,1))
  J <- 2
  T_max <- 4
  apply_lower_bounds <- c(TRUE, FALSE)

  # Call the function
  result <- generate_Ps(Y, J, T_max)

  # Check that the result is a list
  expect_type(result, "list")

  # Check that the length of the list is equal to J
  expect_equal(length(result), J)

  # Check that each element of the list is a matrix
  for (matrices_j in result) {
    expect_type(matrices_j, "list")
    expect_equal(length(matrices_j), T_max-1)
    for (t in 1:(T_max-1)) {
      for (l in 1:length(Y$values)) {
        expect_equal(sum(matrices_j[[t]][l,]), 1)
      }
    }
  }
})


test_that("generate_transition_model returns correct components for treated_period > 0", {
  set.seed(1234)
  J <- 2
  T_max <- 10
  treated_period <- 5
  outcomes_list <- list(1, 1:2, 1:3, 1:4)

  for (outcomes in outcomes_list) {

    result <- generate_transition_model(outcomes, J, T_max, treated_period = treated_period)

    # Check if a valid transition model is returned
    expect_true(is.list(result))
    expect_true("Ps_control" %in% names(result))
    expect_true("Ps_treated" %in% names(result))
    expect_true("pmfs_initial_control" %in% names(result))
    expect_true("priors" %in% names(result))
    expect_true("priors_treated" %in% names(result))
    expect_true("treated_period" %in% names(result))

    # Check the structure of Ps_treated
    expect_true(is.list(result$Ps_treated))
    expect_equal(length(result$Ps_treated), J)

    # Check if first treated_period-1 matrices are identical
    pre_treatment_indices <- 1:(treated_period-2)
    for (t in pre_treatment_indices) {
      expect_true(all(result$Ps_treated[[1]][[t]] ==
                        result$Ps_control[[1]][[t]]))
      expect_true(all(result$Ps_treated[[2]][[t]] ==
                        result$Ps_control[[2]][[t]]))
    }
  }
})

test_that("generate_transition_model handles different treated_periods correctly", {
  outcomes <- c(0, 1)
  J_max <- 3
  T_max <- 10

  for (J in 1:J_max)
    for (treated_period in 1:T_max) {
      result <- generate_transition_model(outcomes, J, T_max, treated_period = treated_period)
      expect_equal(length(result$Ps_treated), J)
      pre_treatment_indices <- 1:(treated_period-2)
      post_treatment_indices <- (treated_period-1):T_max

      for (j in 1:J) {
        if (treated_period > 2) {
          for (t in pre_treatment_indices) {
            expect_true(all(result$Ps_treated[[j]][[t]] ==
                              result$Ps_control[[j]][[t]]))
          }
        }
      }
    }
})
