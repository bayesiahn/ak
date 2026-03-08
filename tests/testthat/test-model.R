test_that("compute_probs_treated_by_type computes correct probabilities", {
  priors_treated <- c(0.6, 0.4)
  priors <- c(0.5, 0.5)
  prob_treated <- 0.8
  expected_probs <- c(0.96, 0.64) # Expected computed probabilities

  computed_probs <- compute_probs_treated_by_type(priors_treated, priors, prob_treated)

  expect_equal(computed_probs, expected_probs)
})

test_that("compute_probs_treated_by_type throws error for invalid probabilities", {
  priors_treated <- c(1.2, 0.8) # These priors would lead to invalid probabilities
  priors <- c(0.5, 0.5)
  prob_treated <- 0.8

  expect_error(compute_probs_treated_by_type(priors_treated, priors, prob_treated))
})

test_that("get_occurence_rates returns correct occurrence rates", {
  p_y1d <- list(matrix(c(0.3, 0.1,
                         0.5, 0.1), byrow = TRUE, ncol = 2),
                matrix(c(0.2, 0.1,
                         0.1, 0.6), byrow = TRUE, ncol = 2))
  priors <- c(0.4, 0.6)
  J <- length(priors)

  occurance_rate <- get_occurence_rates(p_y1d, priors)

  #
  expected_rates <- matrix(c((0.3 + 0.5)*0.4, (0.2 + 0.1)*0.6,
                             (0.1 + 0.1)*0.4, (0.1 + 0.6)*0.6), byrow = TRUE, ncol = 2,
                           dimnames = list(c("Control", "Treated"), 1:J))

  expect_equal(occurance_rate, expected_rates)
  expect_equal(sum(occurance_rate), 1)
})

test_that("get_expected_outcomes_df returns expected structure and content", {
  transition_model <- generate_transition_model(c(0,1))
  outcomes_df <- get_expected_outcomes_df(transition_model)

  expect_s3_class(outcomes_df, "data.frame")
  expect_true(all(c("j", "t", "y", "g", "occurance_rate") %in% names(outcomes_df)))

  # Verify the dimensions of the dataframe
  expected_rows <- (1+length(transition_model$Ps_control[[1]])) * length(transition_model$priors) * 2 # for both treated and control
  expect_equal(nrow(outcomes_df), expected_rows)
  expect_equal(ncol(outcomes_df), 5)

  # Check if expected outcomes are present in the data frame
  expect_equal(outcomes_df %>% dplyr::filter(t == 1) %>% dplyr::filter(g == 0) %>% dplyr::pull(y),
               c(transition_model$pmfs_initial_control[[1]][2],
                 transition_model$pmfs_initial_control[[2]][2]),
               ignore_attr = TRUE)

  expect_equal(outcomes_df %>% dplyr::filter(t == 1) %>% dplyr::filter(g != 0) %>% dplyr::pull(y),
               c(transition_model$pmfs_initial_treated[[1]][2],
                 transition_model$pmfs_initial_treated[[2]][2]),
               ignore_attr = TRUE)
})

test_that("get_expected_outcomes_df aggregates correctly", {
  transition_model <- generate_transition_model()
  plot_data_df <- get_expected_outcomes_df(transition_model, aggregate = TRUE)

  # With aggregation, check for a single 'j' value indicating aggregation
  expect_true(all(plot_data_df$j == "Aggregate (All observations)"))
  # Verify the reduction in the number of rows due to aggregation
  expect_true(nrow(plot_data_df) < nrow(get_expected_outcomes_df(transition_model)))
})

test_that("get_transition_probability_df returns correct structure", {
  transition_model <- generate_transition_model()
  transition_probs_df <- get_transition_probability_df(transition_model)

  expect_true("y_past" %in% names(transition_probs_df))
  expect_true("y_current" %in% names(transition_probs_df))
  expect_true("probability" %in% names(transition_probs_df))
  expect_true("g" %in% names(transition_probs_df))
  expect_true("j" %in% names(transition_probs_df))
})

test_that("get_transition_probability_df calculates probabilities correctly", {
  transition_model <- generate_transition_model(c("Employed", "Unemployed"))
  transition_probs_df <- get_transition_probability_df(transition_model)

  # Example check for specific transition probability
  # Assuming 'Employed' to 'Employed' transition under control conditions
  control_employed_to_employed <- transition_probs_df %>%
    dplyr::filter(y_past == "Employed", y_current == "Employed", g == 0)

  expect_equal(control_employed_to_employed %>% dplyr::filter(j == 1) %>% dplyr::pull(probability),
               sapply(transition_model$Ps_control[[1]], function(x) x[1,1]))
  expect_equal(control_employed_to_employed %>% dplyr::filter(j == 2) %>% dplyr::pull(probability),
               sapply(transition_model$Ps_control[[2]], function(x) x[1,1]))

  # Assuming 'Employed' to 'Employed' transition under treated conditions
  treated_employed_to_employed <- transition_probs_df %>%
    dplyr::filter(y_past == "Employed", y_current == "Employed", g > 0)

  expect_equal(treated_employed_to_employed %>% dplyr::filter(j == 1) %>% dplyr::pull(probability),
               sapply(transition_model$Ps_treated[[1]], function(x) x[1,1]))
  expect_equal(treated_employed_to_employed %>% dplyr::filter(j == 2) %>% dplyr::pull(probability),
               sapply(transition_model$Ps_treated[[2]], function(x) x[1,1]))
})

test_that("modify_outcomes_of_interest correctly updates model outcomes", {
  model <- list(Y = list(names = c("outcome1", "outcome2", "outcome3")))
  model$bootstrap_estimates <- list(model, model)

  # Test case 1: Single valid outcome
  new_model <- modify_outcomes_of_interest(model, c("outcome1"))
  expect_equal(new_model$Y$names, model$Y$names)
  expect_equal(new_model$Y$values, c(1,0,0))
  for (bootstrap_estimate in new_model$bootstrap_estimates) {
    expect_equal(bootstrap_estimate$Y$values, new_model$Y$values)
  }

  # Test case 2: Multiple valid outcomes
  new_model <- modify_outcomes_of_interest(model, c("outcome1", "outcome3"))
  expect_equal(new_model$Y$names, model$Y$names)
  expect_equal(new_model$Y$values, c(1,0,1))
  for (bootstrap_estimate in new_model$bootstrap_estimates) {
    expect_equal(bootstrap_estimate$Y$values, new_model$Y$values)
  }

  # Test case 3: No matching outcomes
  new_model <- modify_outcomes_of_interest(model, c("outcomeX"))
  expect_equal(new_model$Y$names, model$Y$names)
  expect_equal(new_model$Y$values, c(0,0,0))
  for (bootstrap_estimate in new_model$bootstrap_estimates) {
    expect_equal(bootstrap_estimate$Y$values, new_model$Y$values)
  }

  # Test case 4: All outcomes match
  new_model <- modify_outcomes_of_interest(model, c("outcome1", "outcome2", "outcome3"))
  expect_equal(new_model$Y$names, model$Y$names)
  expect_equal(new_model$Y$values, c(1,1,1))
  for (bootstrap_estimate in new_model$bootstrap_estimates) {
    expect_equal(bootstrap_estimate$Y$values, new_model$Y$values)
  }
})


