test_that("get_weighted_transition_probabilities_by_gjt works", {
  set.seed(1234)
  N <- 5000
  g <- sample(c(0, 1), size = 5000, replace = TRUE)
  y <- matrix(sample(c(0, 1), size = length(g)*2, replace = TRUE), nrow = length(g))
  weights <- t(sapply(1:N, function(i) sample(c(0, 1), size = 2)))
  weighted_probs_df <- get_weighted_transition_probabilities_by_gjt(y, g, weights, 1, 1)

  # For the first periods, weighted probabilities are NaN since past values do not exist
  expect_true(all(is.nan(weighted_probs_df %>% dplyr::filter(t == 1) %>% dplyr::pull(weighted_mean))))

  # For the other periods, weighted probabilities should converge around 0.5
  # This test checks if the calculated weighted probabilities for periods other than the first one
  # are approximately equal to 0.5 (within a tolerance of 0.05). This is expected under the assumption
  # that the random sample of y values and weights should average out to around 0.5 over time.
  weighted_probs <- weighted_probs_df %>% dplyr::filter(t != 1) %>% dplyr::pull(weighted_mean)
  expect_equal(weighted_probs, rep(0.5, length(weighted_probs)), tolerance = 0.05)
})

test_that("plot_counterfactuals_comparison_differences returns ggplot object", {
  set.seed(1234)
  # Generate a simple model for testing
  transition_model <- generate_transition_model(c(1, 2), J = 2, T = 4, treated_period = 3)
  sample <- generate_sample(transition_model, 100)
  y <- sample$y
  g <- sample$g

  # Estimate model
  model <- estimate_transition_model_from_wide_matrix(
    J = 2, y = y, g = g,
    opts = list(multistart_counts = 2, bootstrap_counts = 0)
  )

  # Test that the function returns without error and produces a ggplot-like object
  p <- plot_counterfactuals_comparison_differences(model)

  # Check that result is a valid ggplot/gtable object (ggarrange returns gtable)
  expect_true(inherits(p, c("gg", "ggplot", "gtable", "ggarrange")))
})

test_that("bootstrap_estimates_to_ci_df produces symmetric CIs", {
  # Create simple estimates and bootstrap samples
  estimate <- c(a = 1.0, b = 2.0, c = 3.0)
  set.seed(42)
  bootstrap_estimates <- lapply(1:100, function(i) {
    estimate + rnorm(3, sd = 0.5)
  })

  # Test pointwise CIs
  ci_pointwise <- bootstrap_estimates_to_ci_df(estimate, bootstrap_estimates,
                                               uniform_ci = FALSE)

  # CIs should be symmetric: (estimate - ci_lower) == (ci_upper - estimate)
  lower_dist <- ci_pointwise$estimate - ci_pointwise$ci_lower
  upper_dist <- ci_pointwise$ci_upper - ci_pointwise$estimate
  expect_equal(lower_dist, upper_dist, tolerance = 1e-10)

  # Test uniform CIs
  ci_uniform <- bootstrap_estimates_to_ci_df(estimate, bootstrap_estimates,
                                             uniform_ci = TRUE)

  # Uniform CIs should also be symmetric
  lower_dist_uniform <- ci_uniform$estimate - ci_uniform$ci_lower
  upper_dist_uniform <- ci_uniform$ci_upper - ci_uniform$estimate
  expect_equal(lower_dist_uniform, upper_dist_uniform, tolerance = 1e-10)

  # Uniform CIs should use the same critical value for all periods
  # So the half-widths relative to SD should be constant
  half_widths_uniform <- upper_dist_uniform / ci_uniform$sd
  expect_true(all(abs(half_widths_uniform - half_widths_uniform[1]) < 1e-10))
})

test_that("uniform CIs are at least as wide as pointwise CIs", {
  # Create estimates and bootstrap samples
  estimate <- c(a = 1.0, b = 2.0, c = 3.0)
  set.seed(42)
  bootstrap_estimates <- lapply(1:100, function(i) {
    estimate + rnorm(3, sd = 0.5)
  })

  ci_pointwise <- bootstrap_estimates_to_ci_df(estimate, bootstrap_estimates,
                                               uniform_ci = FALSE)
  ci_uniform <- bootstrap_estimates_to_ci_df(estimate, bootstrap_estimates,
                                             uniform_ci = TRUE)

  # Uniform critical value should be >= all pointwise critical values
  expect_true(all(ci_uniform$critical_value >= ci_pointwise$critical_value))

  # Uniform CI widths should be >= pointwise CI widths
  width_pointwise <- ci_pointwise$ci_upper - ci_pointwise$ci_lower
  width_uniform <- ci_uniform$ci_upper - ci_uniform$ci_lower
  expect_true(all(width_uniform >= width_pointwise - 1e-10))
})

test_that("attach_ci_to_plot_data_df works with uniform_ci parameter", {
  set.seed(1234)
  # Generate a simple model for testing
  transition_model <- generate_transition_model(c(1, 2), J = 1, T = 5,
                                                 treated_period = 3)
  sample <- generate_sample(transition_model, 200)
  y <- sample$y
  g <- sample$g

  # Estimate model with bootstrap (se_method = "bootstrap" is required)
  model <- estimate_transition_model_from_wide_matrix(
    J = 1, y = y, g = g,
    opts = list(multistart_counts = 2, se_method = "bootstrap",
                bootstrap_counts = 10)
  )

  # Check that bootstrap_estimates exists
  expect_true("bootstrap_estimates" %in% names(model))

  # Get LTATTs from model summary
  model_summary <- summarize_transition_model(model)
  ltatts <- model_summary$LTATTs[[1]]

  # Create a simple plot data frame matching LTATTs
  # LTATTs start from treated_period (3) to T_max (5)
  post_treatment_periods <- 3:5
  plot_data_df <- data.frame(
    j = rep(1, length(post_treatment_periods)),
    g = rep(3, length(post_treatment_periods)),
    t = post_treatment_periods,
    y = ltatts
  )

  # Test with uniform_ci = TRUE
  result_uniform <- attach_ci_to_plot_data_df(
    plot_data_df, model, "LTATTs", uniform_ci = TRUE
  )

  # Test with uniform_ci = FALSE
  result_pointwise <- attach_ci_to_plot_data_df(
    plot_data_df, model, "LTATTs", uniform_ci = FALSE
  )

  # Both should have ci_lower and ci_upper columns
  expect_true("ci_lower" %in% names(result_uniform))
  expect_true("ci_upper" %in% names(result_uniform))
  expect_true("ci_lower" %in% names(result_pointwise))
  expect_true("ci_upper" %in% names(result_pointwise))

  # CIs should be symmetric around the point estimate (y)
  lower_dist_uniform <- result_uniform$y - result_uniform$ci_lower
  upper_dist_uniform <- result_uniform$ci_upper - result_uniform$y
  expect_equal(lower_dist_uniform, upper_dist_uniform, tolerance = 1e-10)

  lower_dist_pointwise <- result_pointwise$y - result_pointwise$ci_lower
  upper_dist_pointwise <- result_pointwise$ci_upper - result_pointwise$y
  expect_equal(lower_dist_pointwise, upper_dist_pointwise, tolerance = 1e-10)
})
