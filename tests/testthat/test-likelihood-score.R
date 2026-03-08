library(testthat)

# Function to compute numerical derivatives
numerical_derivative <- function(f, x) {
  x_length <- length(x)
  h_abs <- 1e-6
  hs <- lapply(1:x_length, function(k) {
    h_vec <- rep(0, x_length)
    h_vec[k] <- h_abs
    h_vec
  })

  sapply(hs, function(h) (f(x + h) - f(x - h)) / (2 * h_abs))
}

test_that("score_pmf_j returns correct dimensions", {
  N <- 5  # Number of observations
  K <- 3  # Number of outcome categories

  set.seed(42)
  y_indices <- sample(1:K, N, replace = TRUE)  # Random outcome categories
  initial_outcome_dist_j <- runif(K)  # Random probabilities
  initial_outcome_dist_j <- initial_outcome_dist_j / sum(initial_outcome_dist_j)  # Normalize to sum 1
  posteriors_j <- runif(N)  # Random posterior probabilities

  for (minimal in c(TRUE, FALSE)) {
    score_j <- score_pmf_j(y_indices, initial_outcome_dist_j, posteriors_j, minimal)
    expect_equal(dim(score_j), c(N, K-minimal))  # Check dimensions
  }
})

test_that("score_pmf_j is close to numerical derivative", {
  N <- 50
  K <- 3

  set.seed(42)
  y_indices <- sample(1:K, N, replace = TRUE)
  initial_outcome_dist_j <- runif(K)
  initial_outcome_dist_j <- initial_outcome_dist_j / sum(initial_outcome_dist_j)
  posteriors_j <- rep(1, N)

  # Function to compute log likelihood (for numDeriv)
  log_likelihood_fn <- function(p) {
    if (any(p < 0) || sum(p) > 1) return(-Inf)  # Ensure valid probability vector
    adjusted_p <- c(p, 1 - sum(p))  # Ensure last probability is constrained
    sum(log(adjusted_p[y_indices]))  # Log likelihood
  }

  # Compute numerical gradient
  score_j_numerical <- numerical_derivative(log_likelihood_fn, initial_outcome_dist_j[1:(K-1)])

  # Compute analytical score function
  score_j <- score_pmf_j(y_indices, initial_outcome_dist_j, posteriors_j)
  score_j <- apply(score_j[, 1:(K-1)], 2, sum)

  # Compare numerical and analytical results
  expect_true(max(abs(score_j - score_j_numerical)) < 1e-4)
})

test_that("score_pmf_j handles single category (K=1) correctly", {
  N <- 10  # Number of observations
  K <- 1   # Only one outcome category

  set.seed(42)
  y_indices <- rep(1, N)  # All outcomes belong to category 1
  initial_outcome_dist_j <- c(1)  # Only one category, probability 1
  posteriors_j <- runif(N)  # Random posterior probabilities


  for (minimal in c(TRUE, FALSE)) {
    score_j <- score_pmf_j(y_indices, initial_outcome_dist_j, posteriors_j, minimal)
    expect_equal(dim(score_j), c(N, K-minimal))  # Check dimensions
    expect_equal(score_j, matrix(0, nrow = N, ncol = (K-minimal))) # Should be all zeros (no variation)
  }

})


test_that("score_pmfs returns correct dimensions", {
  set.seed(123)
  N <- 100
  J <- 2
  K <- 3
  y_indices <- sample(1:K, N, replace = TRUE)
  pmfs <- list(runif(K), runif(K))
  pmfs <- lapply(pmfs, function(p) p / sum(p))  # Normalize to sum to 1
  posteriors <- matrix(runif(N * J), nrow = N, ncol = J)
  posteriors <- posteriors / rowSums(posteriors)  # Normalize rows


  for (minimal in c(TRUE, FALSE)) {
    scores <- score_pmfs(y_indices, pmfs, posteriors, minimal = minimal)

    expect_equal(length(scores), J)  # Should return a list of J elements
    expect_equal(dim(scores[[1]]), c(N, K-minimal))  # Each element should be N x (K-minimal)
  }
})

test_that("score_pmfs handles valid inputs correctly", {
  set.seed(42)
  N <- 50
  J <- 3
  K <- 4
  y_indices <- sample(1:K, N, replace = TRUE)
  pmfs <- list(runif(K), runif(K), runif(K))
  pmfs <- lapply(pmfs, function(p) p / sum(p))
  posteriors <- matrix(runif(N * J), nrow = N, ncol = J)
  posteriors <- posteriors / rowSums(posteriors)


  for (minimal in c(TRUE, FALSE)) {
    scores <- score_pmfs(y_indices, pmfs, posteriors, minimal = minimal)

    expect_true(all(!is.na(unlist(scores))))  # No NA values
    if (!minimal) {
      expect_true(all(sapply(scores, function(s) all(abs(rowSums(s) - 1) < 1e-6))))  # Ensure sum constraint holds
    }
  }
})

test_that("score_pmfs handles J=1 case correctly", {
  set.seed(555)
  N <- 10
  J <- 1
  K <- 3
  y_indices <- sample(1:K, N, replace = TRUE)
  pmfs <- list(runif(K))
  pmfs <- lapply(pmfs, function(p) p / sum(p))
  posteriors <- matrix(1, nrow = N, ncol = J)


  for (minimal in c(TRUE, FALSE)) {
    scores <- score_pmfs(y_indices, pmfs, posteriors, minimal = minimal)

    expect_equal(length(scores), J)  # Should return a single matrix in a list
    expect_equal(dim(scores[[1]]), c(N, K-minimal))  # Should have correct dimensions
  }
})

test_that("score_pmfs handles K=1 case correctly", {
  set.seed(42)
  N <- 10
  J <- 2
  K <- 1  # Single outcome category
  y_indices <- rep(1, N)
  pmfs <- list(c(1), c(1))  # Probability is always 1 for K=1
  posteriors <- matrix(runif(N * J), nrow = N, ncol = J)
  posteriors <- posteriors / rowSums(posteriors)


  for (minimal in c(TRUE, FALSE)) {
    scores <- score_pmfs(y_indices, pmfs, posteriors, minimal = minimal)

    expect_equal(dim(scores[[1]]), c(N, K-minimal))  # Expect dimensions to be N x (1-minimal)
    expect_equal(unname(scores[[1]]), matrix(0, nrow = N, ncol = (K-minimal)))  # Should return all zeros
  }
})

test_that("score_P_jt returns correct dimensions", {
  K <- 2  # Number of outcome categories
  t <- 1  # Time step

  model <- generate_transition_model(outcomes = 1:K, J = 1, T_max = 3)
  sample <- generate_sample(model, N = 100)
  P_jt <- model$Ps_control[[1]][[t]]  # Transition probabilities
  y_indices_matrix <- sample$y[,1:2]  # Observed outcomes
  posteriors_j <- rep(1, nrow(y_indices_matrix))

  for (minimal in c(TRUE, FALSE)) {
    scores <- score_P_jt(y_indices_matrix, P_jt, posteriors_j, minimal = minimal)

    expect_equal(length(scores), K)
    expect_equal(dim(scores[[1]]), c(nrow(y_indices_matrix), K-minimal))  # Check dimensions
  }
})

test_that("score_P_jt handles valid inputs correctly", {
  set.seed(12)
  K <- 3
  t <- 1

  model <- generate_transition_model(outcomes = 1:K, J = 1, T_max = 3)
  sample <- generate_sample(model, N = 50)
  P_jt <- model$Ps_control[[1]][[t]]
  y_indices_matrix <- sample$y[,1:2]
  posteriors_j <- runif(nrow(y_indices_matrix))


  for (minimal in c(TRUE, FALSE)) {
    scores <- score_P_jt(y_indices_matrix, P_jt, posteriors_j, minimal = minimal)
    scores <- do.call(cbind, scores) # Make it a matrix for easier comparison

    expect_true(all(!is.na(scores)))  # No NA values

    # Check if each row has at least a non-zero value
    expect_true(all(apply(scores, 1, function(row) any(row != 0))))
  }
})


test_that("score_P_jt handles one-single input correctly", {
  set.seed(123)
  K <- 3
  t <- 1
  N <- 1

  model <- generate_transition_model(outcomes = 1:K, J = 1, T_max = 3)
  sample <- generate_sample(model, N = N)
  P_jt <- model$Ps_control[[1]][[t]]
  y_indices_matrix <- matrix(sample$y[,1:2], nrow = N)
  posteriors_j <- runif(nrow(y_indices_matrix))


  for (minimal in c(TRUE, FALSE)) {
    scores <- score_P_jt(y_indices_matrix, P_jt, posteriors_j, minimal = minimal)
    scores <- do.call(cbind, scores) # Make it a matrix for easier comparison

    expect_true(all(!is.na(scores)))  # No NA values

    # Check if each row has at least a non-zero value
    expect_true(all(apply(scores, 1, function(row) any(row != 0))))
  }
})

test_that("score_P_jt is close to numerical derivative", {
  K <- 2
  t <- 1

  model <- generate_transition_model(outcomes = 1:K, J = 1, T_max = 3)
  sample <- generate_sample(model, N = 50)
  P_jt <- matrix(c(0.25, 0.85, 0.75, 0.15), ncol = 2)
  y_indices_matrix <- sample$y[,1:2]
  posteriors_j <- rep(1, nrow(y_indices_matrix))

  # Log-likelihood function
  log_likelihood_fn <- function(P_vec) {
    P_matrix <- matrix(P_vec, nrow = K*(K-1))
    P_matrix <- cbind(P_matrix, 1 - rowSums(P_matrix))  # Ensure each row sums up to 1
    log_likelihoods <- apply(y_indices_matrix, 1, function(y) {
      P_matrix[y[1], y[2]]
    })
    sum(log(log_likelihoods))
  }

  # Compute numerical gradient
  P_vec <- as.vector(P_jt[,1:(K-1)])
  numerical_score <- numerical_derivative(log_likelihood_fn, P_vec)
  numerical_score

  # Compute analytical score function
  score_j <- score_P_jt(y_indices_matrix, P_jt, posteriors_j, minimal = FALSE)
  score_j <- do.call(cbind, score_j) # Make it a matrix for easier comparison
  valid_cols <- 1:(K*K)
  valid_cols <- valid_cols[valid_cols %% K != 0]
  analytical_score <- apply(score_j[,valid_cols], 2, sum)

  # Compare numerical and analytical results
  expect_true(max(abs(analytical_score - numerical_score)) < 1e-4)
})

test_that("score_P_jt handles edge cases (K=1 and probabilities at 0 or 1)", {
  K <- 1
  t <- 1

  model <- generate_transition_model(outcomes = 1:K, J = 1, T_max = 3)
  sample <- generate_sample(model, N = 10)

  for (minimal in c(TRUE, FALSE)) {
    P_jt <- matrix(1, nrow = K, ncol = K)  # Probabilities are all 1
    y_indices_matrix <- sample$y[,1:2]
    posteriors_j <- rep(1, nrow(y_indices_matrix))

    scores <- score_P_jt(y_indices_matrix, P_jt, posteriors_j, minimal = minimal)
    scores <- do.call(cbind, scores) # Make it a matrix for easier comparison
    if (!minimal) {
      expect_equal(scores, matrix(0, nrow = nrow(y_indices_matrix),
                                  ncol = 1, dimnames = list(c(), "1_1")))  # Expect all zeros
    } else {
      expect_equal(scores, matrix(0, nrow = nrow(y_indices_matrix),
                                  ncol = 0))  # Expect empty vector
    }


    # Edge case: One probability is exactly 0
    P_jt <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)  # 2x2 transition matrix
    sample <- generate_sample(model, N = 10)
    y_indices_matrix <- sample$y[,1:2]

    scores <- score_P_jt(y_indices_matrix, P_jt, posteriors_j, minimal = minimal)
    scores <- do.call(cbind, scores) # Make it a matrix for easier comparison
    expect_false(any(is.nan(scores)))  # No NaNs should appear
  }
})


test_that("score_P_j returns correct dimensions", {
  outcomes <- 1:2
  T_max <- 3
  model <- generate_transition_model(outcomes = outcomes, J = 1, T_max = T_max)
  sample <- generate_sample(model, N = 100)
  P_j <- model$Ps_control[[1]]  # List of transition matrices
  y_indices_matrix <- sample$y  # Full sequence of outcomes
  posteriors_j <- rep(1, nrow(y_indices_matrix))


  for (minimal in c(TRUE, FALSE)) {
    scores <- score_P_j(y_indices_matrix, P_j, posteriors_j, minimal = minimal)

    expect_equal(length(scores), (T_max-1))
    expect_equal(length(scores[[1]]), length(outcomes))
    expect_equal(dim(scores[[1]][[1]]), c(nrow(y_indices_matrix), length(outcomes)-minimal))
  }
})

test_that("score_P_j handles valid inputs correctly", {
  set.seed(42)
  outcomes <- 1:3
  T_max <- 2
  model <- generate_transition_model(outcomes = outcomes, J = 1, T_max = T_max)
  sample <- generate_sample(model, N = 50)
  P_j <- model$Ps_control[[1]]
  y_indices_matrix <- sample$y
  posteriors_j <- runif(nrow(y_indices_matrix))

  for (minimal in c(TRUE, FALSE)) {
    scores <- score_P_j(y_indices_matrix, P_j, posteriors_j, minimal = minimal)
    scores <- lapply(scores, function(scores_l) do.call(cbind, scores_l)) # Make it a matrix for easier comparison
    scores <- do.call(cbind, scores)

    expect_true(all(!is.na(scores)))  # No NA values
    expect_true(all(apply(scores, 1, function(row) any(row != 0))))  # Ensure non-zero values
  }
})

test_that("score_P_j is close to numerical derivative", {
  set.seed(555)
  K <- 2
  outcomes <- 1:K
  T_max <- 3
  model <- generate_transition_model(outcomes = outcomes, J = 1, T_max = T_max)
  sample <- generate_sample(model, N = 5)
  P_j <- list(matrix(c(0.4, 0.7, 0.6, 0.3), ncol = 2), matrix(c(0.5, 0.8, 0.5, 0.2), ncol = 2))
  y_indices_matrix <- sample$y
  posteriors_j <- rep(1, nrow(y_indices_matrix))

  # Log-likelihood function
  log_likelihood_fn <- function(P_vec) {
    P_matrix_1 <- matrix(P_vec[1:(K*(K-1))], nrow = K*(K-1))
    P_matrix_2 <- matrix(P_vec[(K*(K-1)+1):length(P_vec)], nrow = K*(K-1))
    P_matrix_1 <- cbind(P_matrix_1, 1 - rowSums(P_matrix_1))  # Ensure each row sums up to 1
    P_matrix_2 <- cbind(P_matrix_2, 1 - rowSums(P_matrix_2))  # Ensure each row sums up to 1
    P_list <- list(P_matrix_1, P_matrix_2)
    log_likelihoods <- sapply(1:nrow(y_indices_matrix), function(i) {
      sum(log(sapply(1:(T_max-1), function(t) P_list[[t]][y_indices_matrix[i, t], y_indices_matrix[i, t+1]])))
    })
    sum((log_likelihoods))
  }

  # Compute numerical gradient
  P_vec <- c(as.vector(P_j[[1]][,1:(K-1)]), as.vector(P_j[[2]][,1:(K-1)]))
  numerical_score <- numerical_derivative(log_likelihood_fn, P_vec)
  numerical_score

  # Compute analytical score function
  score_j <- score_P_j(y_indices_matrix, P_j, posteriors_j, minimal = FALSE)
  score_j <- lapply(score_j, function(scores_l) do.call(cbind, scores_l)) # Make it a matrix for easier comparison
  score_j <- do.call(cbind, score_j)

  analytical_score <- colSums(score_j)
  valid_cols <- 1:(K*K*(T_max-1))
  valid_cols <- valid_cols[valid_cols %% K != 0]
  analytical_score <- apply(score_j[,valid_cols], 2, sum)
  analytical_score

  # Compare numerical and analytical results
  expect_true(max(abs(analytical_score - numerical_score)) < 1e-4)
})

test_that("score_P_j handles edge cases (T=1 and probabilities at 0 or 1)", {
  outcomes <- 1:2
  T_max <- 2
  model <- generate_transition_model(outcomes = outcomes, J = 1, T_max = T_max)
  sample <- generate_sample(model, N = 10)
  y_indices_matrix <- sample$y
  posteriors_j <- rep(1, nrow(y_indices_matrix))

  # Edge case: One probability is exactly 0
  P_j <- list(matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE))  # 2x2 transition matrix
  sample <- generate_sample(model, N = 10)
  y_indices_matrix <- sample$y

  for (minimal in c(TRUE, FALSE)) {
    scores <- score_P_j(y_indices_matrix, P_j, posteriors_j, minimal = minimal)
    scores <- lapply(scores, function(scores_l) do.call(cbind, scores_l)) # Make it a matrix for easier comparison
    scores <- do.call(cbind, scores)
    expect_false(any(is.nan(scores)))  # No NaNs should appear
  }
})

test_that("score_Ps returns correct dimensions", {
  outcomes <- 1:2
  T_max <- 3
  J <- 2
  model <- generate_transition_model(outcomes = outcomes, J = J, T_max = T_max)
  sample <- generate_sample(model, N = 100)
  Ps <- model$Ps_control  # List of transition matrices for each component
  y_indices_matrix <- sample$y  # Observed outcomes
  posteriors <- transitions_to_posteriors(sample$y,model$pmfs_initial_control,
                                          model$Ps_control, model$priors)

  for (minimal in c(TRUE, FALSE)) {
    scores <- score_Ps(y_indices_matrix, Ps, posteriors, minimal = minimal)

    expect_equal(length(scores), (J))
    expect_equal(length(scores[[1]]), (T_max-1))
    expect_equal(length(scores[[1]][[1]]), length(outcomes))
    expect_equal(dim(scores[[1]][[1]][[1]]), c(nrow(y_indices_matrix), length(outcomes)-minimal))
  }
})

test_that("score_Ps handles valid inputs correctly", {
  set.seed(123)
  outcomes <- 1:2
  T_max <- 2
  J <- 2
  model <- generate_transition_model(outcomes = outcomes, J = J, T_max = T_max)
  sample <- generate_sample(model, N = 50)
  Ps <- model$Ps_control
  y_indices_matrix <- sample$y
  posteriors <- transitions_to_posteriors(sample$y,model$pmfs_initial_control,
                                          model$Ps_control, model$priors)

  for (minimal in c(TRUE, FALSE)) {
    scores <- score_Ps(y_indices_matrix, Ps, posteriors, minimal = minimal)
    scores <- lapply(scores, function(scores_j) do.call(cbind,
      lapply(scores_j, function(scores_jt) do.call(cbind, scores_jt))))
    scores <- do.call(cbind, scores) # Make it a matrix for easier comparison

    expect_true(all(!is.na(scores)))  # No NA values
    expect_true(all(apply(scores, 1, function(row) any(row != 0))))  # Ensure non-zero values
  }
})

test_that("score_Ps is close to numerical derivative", {
  set.seed(555)
  K <- 2
  T_max <- 3
  J <- 2
  outcomes <- 1:K
  model <- generate_transition_model(outcomes = outcomes, J = J, T_max = T_max)
  sample <- generate_sample(model, N = 5)
  N <- nrow(sample$y)
  Ps <- list(
    list(matrix(c(0.4, 0.7, 0.6, 0.3), ncol = 2), matrix(c(0.5, 0.8, 0.5, 0.2), ncol = 2)),
    list(matrix(c(0.2, 0.6, 0.8, 0.4), ncol = 2), matrix(c(0.6, 0.3, 0.4, 0.7), ncol = 2))
  )
  y_indices_matrix <- sample$y
  posteriors <- matrix(1, nrow = N, ncol = J)

  # Log-likelihood function
  log_likelihood_fn <- function(P_vec) {
    P_list_1 <- list(matrix(P_vec[1:(K*(K-1))], nrow = K*(K-1)),
                     matrix(P_vec[(K*(K-1)+1):(2*K*(K-1))], nrow = K*(K-1)))
    P_list_2 <- list(matrix(P_vec[(2*K*(K-1)+1):(3*K*(K-1))], nrow = K*(K-1)),
                     matrix(P_vec[(3*K*(K-1)+1):length(P_vec)], nrow = K*(K-1)))

    P_list_1 <- lapply(P_list_1, function(P) cbind(P, 1 - rowSums(P)))
    P_list_2 <- lapply(P_list_2, function(P) cbind(P, 1 - rowSums(P)))

    Ps_test <- list(P_list_1, P_list_2)
    log_likelihoods <- sapply(1:nrow(y_indices_matrix), function(i) {
      sum((sapply(1:J, function(j) {
        sum(log(sapply(1:(T_max-1), function(t) Ps_test[[j]][[t]][y_indices_matrix[i, t], y_indices_matrix[i, t+1]])))
      })))
    })
    sum(log_likelihoods)
  }

  # Compute numerical gradient
  P_vec <- c(as.vector(Ps[[1]][[1]][,1:(K-1)]), as.vector(Ps[[1]][[2]][,1:(K-1)]),
             as.vector(Ps[[2]][[1]][,1:(K-1)]), as.vector(Ps[[2]][[2]][,1:(K-1)]))
  numerical_score <- numerical_derivative(log_likelihood_fn, P_vec)
  numerical_score

  # Compute analytical score function
  score_j <- score_Ps(y_indices_matrix, Ps, posteriors, minimal = FALSE)
  score_j <- lapply(score_j, function(scores_j) do.call(cbind,
                                                      lapply(scores_j, function(scores_jt) do.call(cbind, scores_jt))))
  score_j <- do.call(cbind, score_j) # Make it a matrix for easier comparison
  analytical_score <- colSums(score_j)
  valid_cols <- 1:ncol(score_j)
  valid_cols <- valid_cols[valid_cols %% K != 0]
  analytical_score <- apply(score_j[,valid_cols], 2, sum)
  analytical_score

  # Compare numerical and analytical results
  expect_true(max(abs(analytical_score - numerical_score)) < 1e-6)
})

test_that("score_priors returns correct dimensions", {
  N <- 100
  J <- 3
  priors <- c(0.3, 0.5, 0.2)
  posteriors <- matrix(runif(N * J), nrow = N, ncol = J)
  posteriors <- posteriors / rowSums(posteriors)  # Normalize rows to sum to 1

  for (minimal in c(TRUE, FALSE)) {
    scores <- score_priors(priors, posteriors, minimal)

    expect_equal(dim(scores), c(N, J-minimal))  # Check dimensions
  }
})

test_that("score_priors handles valid inputs correctly", {
  set.seed(42)
  N <- 50
  J <- 2
  priors <- c(0.6, 0.4)
  posteriors <- matrix(runif(N * J), nrow = N, ncol = J)
  posteriors <- posteriors / rowSums(posteriors)  # Normalize rows

  for (minimal in c(TRUE, FALSE)) {
    scores <- score_priors(priors, posteriors, minimal = minimal)

    expect_true(all(!is.na(scores)))  # No NA values
  }
})

test_that("score_priors handles J=1 case correctly", {
  N <- 10
  J <- 1
  priors <- c(1)  # Only one component, must have probability 1
  posteriors <- matrix(1, nrow = N, ncol = J)  # Posterior is always 1 for J=1

  for (minimal in c(TRUE, FALSE)) {
    scores <- score_priors(priors, posteriors, minimal = minimal)

    expect_equal(unname(scores), matrix(0, nrow = N, ncol = (J-minimal)))  # Should be all zeros
  }
})

test_that("score_priors handles edge cases with probabilities 0 or 1", {
  N <- 10
  J <- 3
  priors <- c(0.4, 0.6, 0)  # One prior is zero
  posteriors <- matrix(runif(N * J), nrow = N, ncol = J)
  posteriors <- posteriors / rowSums(posteriors)  # Normalize rows

  expect_error(score_priors(priors, posteriors),
               "Prior probabilities must be non-zero", fixed = TRUE)

  priors <- c(0.5, 0.5, 1)  # Priors are not summing up to one
  for (minimal in c(TRUE, FALSE)) {
    scores <- score_priors(priors, posteriors, minimal = minimal)

    expect_false(any(is.nan(scores)))  # No NaNs should appear
  }
})
