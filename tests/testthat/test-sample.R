test_that("generate_sample_j returns correct output", {
  # Define inputs
  J_max <- 3
  N_j <- 50
  T_max <- 10
  outcomes_list <- list(c(0,1), c(0,1,2))

  # Check for each possible J and outcomes
  for (outcomes in outcomes_list) {
    for (J in 1:J_max) {
      transition_model <- generate_transition_model(outcomes, J, T_max)
      for (j in 1:J) {
        Ps_control_j <- transition_model$Ps_control[[j]]
        pmfs_initial_control_j <- transition_model$pmfs_initial_control[[j]]
        y_j <- generate_sample_j(transition_model$Y, N_j,
                                 Ps_control_j, pmfs_initial_control_j)

        # Check that the result is a list
        expect_type(y_j, "double")

        # Check that the dimension of y_j
        expect_equal(dim(y_j), c(N_j, T_max))
      }
    }
  }
})

test_that("Large sample checks for generate_sample_j", {
  # Define inputs
  set.seed(1234)
  J_max <- 2
  N_j <- 1000
  T_max <- 10

  # Check for each possible J
  for (J in 1:J_max) {
    transition_model <- generate_transition_model(c(0,1), J, T_max,
                                      treated_period = 0)
    for (j in 1:J) {
      Ps_control_j <- transition_model$Ps_control[[j]]
      pmfs_initial_control_j <- as.numeric(transition_model$pmfs_initial_control[[j]])
      y_j <- generate_sample_j(transition_model$Y, N_j,
                               Ps_control_j, pmfs_initial_control_j)

      # Check that the sampled probabilities are close to the true probabilities
      y_j0 <- table(y_j[,1])
      pmfs_initial_control_sampled <- as.numeric(y_j0 / sum(y_j0))
      expect_equal(pmfs_initial_control_sampled, pmfs_initial_control_j,
                   tolerance = 1e-1)

      # Check that the dimension of y_j
      y_j1 <- table(y_j[,2])
      y_j1_dist_sampled <- as.numeric(y_j1 / sum(y_j1))
      y_j1_dist_true <- as.numeric(pmfs_initial_control_j %*% Ps_control_j[[1]])
      expect_equal(y_j1_dist_sampled, y_j1_dist_true,
                   tolerance = 1e-1)

    }
  }
})


test_that("Large sample checks for generate_sample", {
  # Define inputs
  set.seed(1234)
  J <- 1
  N <- 1000
  T_max <- 3

  # Construct transition_model and sample
  transition_model <- generate_transition_model(c(0,1), J, T_max, treated_period = 0)
  sample <- generate_sample(transition_model, N = N)
  weights <- matrix(1, nrow = N, ncol = 1)
  y_values <- transition_model$Y$values

  for (y_past_index in 1:length(y_values))
    for (y_current_index in 1:length(y_values)) {
      y_past <- y_values[y_past_index]
      y_current <- y_values[y_current_index]

      # Construct sampled probabilities
      df <- get_weighted_transition_probabilities_by_gjt(sample$y, sample$g, weights, y_current, y_past)
      weighted_means <- df %>% dplyr::filter(t != 1) %>% dplyr::pull(weighted_mean)
      true_means <- sapply(transition_model$Ps_control[[1]], function(x) x[y_past_index, y_current_index])

      # Check that the sampled probabilities are close to the true probabilities
      expect_equal(weighted_means, true_means,
                   tolerance = 1e-1)
    }
})



test_that("Large sample checks for generate_sample, J = 2", {
  # Define inputs
  set.seed(1234)
  J <- 2
  N <- 3000
  T_max <- 3

  # Construct transition_model and sample
  transition_model <- generate_transition_model(c(0,1), J = 2, T_max, treated_period = 0)
  sample <- generate_sample(transition_model, N = N)
  y_values <- transition_model$Y$values

  # Construct weights
  weights <- matrix(0, ncol = J, nrow = N)
  weights[sample$j == 1, 1] <- 1
  weights[sample$j == 2, 2] <- 1

  for (y_past_index in 1:length(y_values))
    for (y_current_index in 1:length(y_values)) {
      y_past <- y_values[y_past_index]
      y_current <- y_values[y_current_index]

      # Construct sampled probabilities
      df <- get_weighted_transition_probabilities_by_gjt(sample$y, sample$g, weights, y_current, y_past)
      for (j_select in 1:J) {
        weighted_means <- df %>% dplyr::filter(t != 1) %>%  dplyr::filter(j == j_select) %>% dplyr::pull(weighted_mean)
        true_means <- sapply(transition_model$Ps_control[[j_select]], function(x) x[y_past_index, y_current_index])

        # Check that the sampled probabilities are close to the true probabilities
        expect_equal(weighted_means, true_means,
                     tolerance = 1e-1)
      }
    }
})

test_that("Large sample checks for priors_treated", {
  # Define inputs
  prob_treated_list <- list(0.2, 0.3, 0.4)
  priors_list <- list(c(0.7, 0.3), c(0.5, 0.5))
  priors_treated_list <- list(c(0.3, 0.7), c(0.4, 0.6), c(0.5, 0.5))

  # Define inputs
  set.seed(123456)
  J <- 2
  N <- 3000
  T_max <- 2

  for (prob_treated in prob_treated_list) {
    for (priors in priors_list) {
      for (priors_treated in priors_treated_list) {
        # Construct transition_model and sample
        transition_model <- generate_transition_model(c(0,1), J = 2, T_max)
        transition_model$prob_treated <- prob_treated
        transition_model$priors <- priors
        transition_model$priors_treated <- priors_treated
        sample <- generate_sample(transition_model, N = N)

        # Construct weights
        weights <- matrix(0, ncol = J, nrow = N)
        weights[sample$j == 1, 1] <- 1
        weights[sample$j == 2, 2] <- 1

        # Check priors_treated
        priors_treated_sample <- c(apply(weights[sample$g > 0, ], 2, mean))
        priors_treated_true <- as.numeric(unname(transition_model$priors_treated))

        # Check that the sampled probabilities are close to the true probabilities
        expect_equal(priors_treated_sample, priors_treated_true,
                     tolerance = 1e-1)
      }
    }
  }
})



test_that("generate_sample returns correct output", {
  # Define inputs
  J_max <- 5
  N <- 50
  T_max <- 10
  outcomes_list <- list(c(0,1), c(0,1,2))

  # Check for each possible J and outcomes
  for (outcomes in outcomes_list) {
    for (J in 1:J_max) {
      transition_model <- generate_transition_model(c(0,1), J, T_max)
      sample <- generate_sample(transition_model, N)

      # Check that the result is a list
      expect_type(sample, "list")

      # Check that the list has the correct names
      expect_equal(names(sample), c("y", "g", "j", "y_its_control", "y_its_treated", "N_js", "N_gjs_table", "transition_model"))

      # Check that y is a matrix
      expect_type(sample$y, "double")

      # Check that g is a vector and has the same length as rows of y
      expect_type(sample$g, "double")
      expect_equal(length(sample$g), nrow(sample$y))

      # Check that j is a vector and has the same length as rows of y
      expect_type(sample$g, "double")
      expect_equal(length(sample$g), nrow(sample$y))

      # Check that y_its is a list
      expect_type(sample$y_its_control, "list")
      expect_type(sample$y_its_treated, "list")

      # Check that N_js is a table
      expect_type(sample$N_js, "integer")

      # Check that transition_model is a list
      expect_type(sample$transition_model, "list")
    }
  }
})

test_that("Expect error for invalid priors_treated", {
  J <- 2
  T_max <- 2
  transition_model <- generate_transition_model(c(0,1), J, T_max)
  transition_model$priors_treated <- c(0.8, 0.2)
  transition_model$prob_treated <- 0.5
  transition_model$priors <- c(0.2, 0.8)

  # This gives an invalid priors_treated because
  # P(D_i = 1| Z_i = 1) is P(Z_i = 1 | D_i = 1) * P(D_i = 1) / P(Z_i = 1)
  # OR priors_treated[j] * prob_treated / priors[j] with j = 1
  # which is 0.8 * 0.5 / 0.2 = 2 > 1
  expect_error(generate_sample(transition_model, N = 10))
})

test_that("Empirical frequencies at T_0 match theoretical pre-treatment PMFs", {
  # 1. Setup: Generate a model and a large sample to ensure convergence
  set.seed(123)
  J <- 2
  T_max <- 10
  treated_period <- 7
  N <- 10000 # Large N to reduce sampling error

  model <- generate_transition_model(outcomes = c(0, 1), J = J,
                                     T_max = T_max,
                                     treated_period = treated_period)

  sample_data <- generate_sample(model, N = N)

  # 2. Identify indices and outcomes just before treatment
  # Column index for T_0 is treated_period - 1
  t0_idx <- model$treated_period - 1
  y_t0 <- sample_data$y[, t0_idx]
  g <- sample_data$g
  j_latent <- sample_data$j # Latent group assignments

  # 3. Define a helper function to calculate empirical PMF
  get_empirical_pmf <- function(outcomes, outcome_levels) {
    counts <- table(factor(outcomes, levels = outcome_levels))
    as.numeric(counts / sum(counts))
  }

  outcome_names <- model$Y$names

  # 4. Iterate through each latent group to check consistency
  for (j in 1:J) {
    # Theoretical PMFs from the model
    expected_control_pmf <- model$pmfs_pre_treatment_control[[j]]
    expected_treated_pmf <- model$pmfs_pre_treatment_treated[[j]]

    # Extract empirical outcomes for units in latent group j
    # Control units (g == 0)
    obs_control <- y_t0[j_latent == j & g == 0]
    # Treated units (g == treated_period)
    obs_treated <- y_t0[j_latent == j & g == model$treated_period]

    # Calculate empirical PMFs
    emp_control_pmf <- get_empirical_pmf(obs_control, outcome_names)
    emp_treated_pmf <- get_empirical_pmf(obs_treated, outcome_names)

    # 5. Assertions: Check if absolute differences are within a reasonable tolerance
    # (Using 0.05 tolerance for N=10000; can be tighter with larger N)
    expect_equal(emp_control_pmf, expected_control_pmf, tolerance = 0.05,
                 label = paste("Control empirical PMF for group", j))

    expect_equal(emp_treated_pmf, expected_treated_pmf, tolerance = 0.05,
                 label = paste("Treated empirical PMF for group", j))
  }
})
