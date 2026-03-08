test_that("score_to_pmf returns correct dimensions", {
  set.seed(123)
  N <- 100
  K <- 3

  for (minimal in c(TRUE, FALSE)) {
    pmf_j_score <- matrix(runif(N * (K-minimal)), nrow = N, ncol = (K-minimal))
    weights <- runif(N)

    pmf <- score_to_pmf(pmf_j_score, rep(1:K)/K, weights, minimal)
    expect_equal(length(pmf), K)
  }
})


test_that("score_to_pmfs returns correct dimensions", {
  set.seed(123)
  N <- 100
  K <- 3
  J <- 2

  for (minimal in c(TRUE, FALSE)) {
    pmfs_score <- lapply(1:J, function(x) matrix(runif(N * (K-minimal)), nrow = N, ncol = (K-minimal)))
    weights <- runif(N)

    pmfs <- list(runif(K)/K, runif(K)/K)
    pmfs <- score_to_pmfs(pmfs_score, pmfs, weights, minimal)

    expect_equal(length(pmfs), J)  # Should return J elements
    expect_equal(length(pmfs[[1]]), K)  # Each PMF should have K elements
  }
})

test_that("score_to_pmfs handles missing weights correctly", {
  set.seed(123)
  N <- 50
  K <- 4
  J <- 3

  for (minimal in c(TRUE, FALSE)) {
    pmfs_score <- lapply(1:J, function(x) matrix(runif(N * (K-minimal)), nrow = N, ncol = (K-minimal)))

    pmfs <- list(runif(K)/K, runif(K)/K, runif(K)/K)
    pmfs <- score_to_pmfs(pmfs_score, pmfs, minimal = minimal)

    expect_equal(length(pmfs), J)
    expect_equal(length(pmfs[[1]]), K)
  }
})

test_that("score_to_pmfs handles incorrect weight length", {
  set.seed(123)
  N <- 100
  K <- 3
  J <- 2


  for (minimal in c(TRUE, FALSE)) {
    pmfs_score <- lapply(1:J, function(x) matrix(runif(N * (K-minimal)), nrow = N, ncol = (K-minimal)))
    weights <- runif(N - 1)  # Incorrect length
    pmfs <- list(runif(K)/K, runif(K)/K)

    expect_error(score_to_pmfs(pmfs_score, pmfs, weights, minimal), "Length of weights must be equal to the number of rows")
  }
})

test_that("score_to_P_jt returns correct dimensions", {
  set.seed(123)
  N <- 100
  K_t <- 3
  K_t_minus_1 <- 4


  for (minimal in c(TRUE, FALSE)) {
    P_jt_score <- lapply(1:K_t_minus_1, function(x) matrix(runif(N * (K_t-minimal)), nrow = N))
    weights <- runif(N)

    P_jt <- matrix(runif(K_t * K_t_minus_1), nrow = K_t_minus_1, ncol = K_t)
    # normalize
    for (k in 1:K_t_minus_1) {
      P_jt[k, ] <- P_jt[k, ] / sum(P_jt[k, ])
    }

    P_jt <- score_to_P_jt(P_jt_score, P_jt, weights, minimal)
    expect_equal(dim(P_jt), c(K_t_minus_1, K_t))
  }
})
