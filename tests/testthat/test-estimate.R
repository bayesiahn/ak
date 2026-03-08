library(testthat)

test_that("estimate_transition_model_from_wide_matrix returns correct components", {
  # Setup mock data for J, y, g, and opts
  set.seed(1234)
  N <- 50
  J <- 2
  T_max <- 5
  treated_period <- 5
  opts <- list(hard_classification = FALSE)
  transition_model <- generate_transition_model(outcomes = c(1, 2), J = J, T_max = T_max)
  transition_model$treated_period <- treated_period
  sample <- generate_sample(transition_model, N = N)
  y <- sample$y
  g <- sample$g

  # Call estimate_transition_model_from_wide_matrix and verify the result
  result <- estimate_transition_model_from_wide_matrix(J, y, g, opts = opts)

  # Check that the result is a list and contains expected components
  expect_type(result, "list")
  expect_true("posteriors" %in% names(result))
  expect_true("weights" %in% names(result))
  expect_true("unit_weights" %in% names(result))
  expect_true("Ps_control_empirical" %in% names(result))
  expect_true("Ps_treated_empirical" %in% names(result))

  # Check if weights are all identical as posteriors, as hard_classification is off
  expect_equal(result$weights, result$posteriors)

  # Check if estimate_transition_model returns exactly identical results
  data <- data.frame(y = c(t(y)),
                     id = rep(1:N, each = (T_max)),
                     t = rep(1:(T_max), N),
                     g = rep(g, each = (T_max)))
  result2 <- estimate_transition_model(J, data, "y", "g", "id", "t", opts = opts)
  expect_equal(result$weights, result2$weights)

  # Pretreatment periods would have slightly different transition probabilities
  for (t in 1:(treated_period-2)) {
    for (j in 1:J) {
      expect_false(all(result$Ps_control[[j]][[t]] == result$Ps_control_empirical[[j]][[t]]))
      expect_false(all(result$Ps_treated[[j]][[t]] == result$Ps_treated_empirical[[j]][[t]]))
    }
  }

  # Post-treatment periods should have identical transition probabilities
  for (t in (treated_period-1):(T_max-1)) {
    for (j in 1:J) {
      expect_equal(result$Ps_control[[j]][[t]], result$Ps_control_empirical[[j]][[t]], tolerance = 1e-2)
      expect_equal(result$Ps_treated[[j]][[t]], result$Ps_treated_empirical[[j]][[t]], tolerance = 1e-2)
    }
  }
})


test_that("estimate_transition_model_from_wide_matrix with hard classification returns correct components", {
  # Setup mock data for J, y, g, and opts
  set.seed(1234)
  J <- 2
  T_max <- 5
  treated_period <- 5
  opts <- list(hard_classification = TRUE)
  transition_model <- generate_transition_model(outcomes = c(1, 2), J = J, T_max = T_max)
  transition_model$treated_period <- treated_period
  sample <- generate_sample(transition_model, N = 50)
  y <- sample$y
  g <- sample$g

  # Call estimate_transition_model_from_wide_matrix and verify the result
  result <- estimate_transition_model_from_wide_matrix(J, y, g, opts = opts)

  # Check that the result is a list and contains expected components
  expect_type(result, "list")
  expect_true("posteriors" %in% names(result))
  expect_true("weights" %in% names(result))
  expect_true("Ps_control_empirical" %in% names(result))
  expect_true("Ps_treated_empirical" %in% names(result))

  # Check if weights are not identical as posteriors, as hard_classification is on
  expect_true(any(c(result$weights != result$posteriors)))

  # Check if rowsums in weights are all 1
  expect_true(all(rowSums(result$weights) == 1))
})

test_that("reorder_transition_model correctly reorders components", {
  # Generate a mock TransitionModel object
  mock_transition_model <- list(
    priors = c(0.3, 0.7),
    priors_treated = c(0.4, 0.6),
    Ps_control = list(matrix(1:4, 2, 2), matrix(5:8, 2, 2)),
    Ps_treated = list(matrix(9:12, 2, 2), matrix(13:16, 2, 2)),
    pmfs_initial_control = list(c(0.1, 0.9), c(0.2, 0.8)),
    pmfs_initial_treated = list(c(0.3, 0.7), c(0.4, 0.6))
  )
  class(mock_transition_model) <- "TransitionModel"

  # Reorder the mock object
  reordered_model <- reorder_transition_model(mock_transition_model)

  # Check if the components are reordered correctly
  expect_equal(reordered_model$priors, sort(mock_transition_model$priors, decreasing = TRUE))
  expect_equal(reordered_model$priors_treated, mock_transition_model$priors_treated[order(mock_transition_model$priors, decreasing = TRUE)])
  expect_equal(reordered_model$Ps_control, mock_transition_model$Ps_control[order(mock_transition_model$priors, decreasing = TRUE)])
  expect_equal(reordered_model$Ps_treated, mock_transition_model$Ps_treated[order(mock_transition_model$priors, decreasing = TRUE)])
  expect_equal(reordered_model$pmfs_initial_control, mock_transition_model$pmfs_initial_control[order(mock_transition_model$priors, decreasing = TRUE)])
  expect_equal(reordered_model$pmfs_initial_treated, mock_transition_model$pmfs_initial_treated[order(mock_transition_model$priors, decreasing = TRUE)])
})

test_that("reorder_transition_model handles non-TransitionModel input", {
  expect_error(reorder_transition_model(list()))
})

# Test 1: Simple matrix
test_that("correctly classifies a simple matrix", {
  mat <- matrix(c(0.1, 0.9, 0.8, 0.2), nrow = 2)
  expect_equal(hard_classification(mat), matrix(c(0, 1, 1, 0), nrow = 2))
  expect_true(all(rowSums(hard_classification(mat)) == 1))
})

# Test 2: Matrix with equal values
test_that("handles ties by selecting the first maximum", {
  mat <- matrix(c(0.5, 0.5, 0.7, 0.7), nrow = 2)
  expect_equal(hard_classification(mat), matrix(c(0, 0, 1, 1), nrow = 2))
  expect_true(all(rowSums(hard_classification(mat)) == 1))
})

# Test 3: Single row matrix
test_that("handles single row matrix", {
  mat <- matrix(c(0.2, 0.3, 0.5), nrow = 1)
  expect_equal(hard_classification(mat), matrix(c(0, 0, 1), nrow = 1))
  expect_true(all(rowSums(hard_classification(mat)) == 1))
})

# Test 4: Single column matrix
test_that("handles single column matrix", {
  mat <- matrix(c(0.2, 0.8), ncol = 1)
  expect_equal(hard_classification(mat), matrix(c(1, 1), ncol = 1))
  expect_true(all(rowSums(hard_classification(mat)) == 1))
})

# Test 5: Large matrix
test_that("handles large matrix", {
  mat <- matrix(runif(100), nrow = 10)
  result <- hard_classification(mat)
  expect_equal(dim(result), c(10, 10))
  expect_true(all(rowSums(result) == 1))
})
