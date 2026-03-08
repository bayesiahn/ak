# tests/testthat/test-summarize_Ps.R
# -------------------------------------------------------------------
# Unit tests for summarize_Ps()
# -------------------------------------------------------------------

test_that("summarize_Ps: basic correctness with defaults", {
  # J = 1, (T-1) = 2, K = 2
  Ps <- list(list(
    matrix(c(0.9, 0.1,
             0.2, 0.8), nrow = 2, byrow = TRUE),
    matrix(c(0.7, 0.3,
             0.4, 0.6), nrow = 2, byrow = TRUE)
  ))

  out <- summarize_Ps(Ps)

  # Expect shape J*(T-1)*K*K = 1*2*2*2 = 8 rows, 5 cols
  expect_equal(nrow(out), 8L)
  expect_equal(ncol(out), 5L)

  # Column names and types
  expect_setequal(names(out), c("j", "t", "from", "to", "prob"))
  expect_type(out$j, "integer")
  expect_type(out$t, "integer")
  expect_type(out$prob, "double")
  expect_true(is.character(out$from))
  expect_true(is.character(out$to))

  # Default y_names should be "0","1" for K=2
  expect_setequal(unique(out$from), c("0", "1"))
  expect_setequal(unique(out$to),   c("0", "1"))

  # Check ordering induced by loops: for a given t,
  # from = 1..K, and for each from, to = 1..K
  # t = 1 block should correspond to first matrix
  expected_t1 <- data.frame(
    j    = 1L,
    t    = 1L,
    from = c("0", "0", "1", "1"),
    to   = c("0", "1", "0", "1"),
    prob = c(0.9, 0.1, 0.2, 0.8),
    stringsAsFactors = FALSE
  )
  expect_equal(out[1:4, ], expected_t1, ignore_attr = TRUE)

  # t = 2 block should correspond to second matrix
  expected_t2 <- data.frame(
    j    = 1L,
    t    = 2L,
    from = c("0", "0", "1", "1"),
    to   = c("0", "1", "0", "1"),
    prob = c(0.7, 0.3, 0.4, 0.6),
    stringsAsFactors = FALSE
  )
  expect_equal(out[5:8, ], expected_t2, ignore_attr = TRUE)
})

test_that("summarize_Ps: respects provided y_names", {
  Ps <- list(list(
    matrix(c(0.5, 0.5,
             0.25, 0.75), nrow = 2, byrow = TRUE)
  ))
  out <- summarize_Ps(Ps, y_names = c("A", "B"))
  expect_setequal(unique(out$from), c("A", "B"))
  expect_setequal(unique(out$to),   c("A", "B"))

  # Spot check probabilities preserved
  # from A -> to B is entry (1,2) = 0.5
  expect_equal(
    out$prob[out$from == "A" & out$to == "B" & out$t == 1 & out$j == 1],
    0.5
  )
  # from B -> to A is entry (2,1) = 0.25
  expect_equal(
    out$prob[out$from == "B" & out$to == "A" & out$t == 1 & out$j == 1],
    0.25
  )
})

test_that("summarize_Ps: handles multiple classes (J) and periods (T-1)", {
  # J = 2, (T-1) = 2, K = 3
  P1_t1 <- matrix(c(
    0.8, 0.1, 0.1,
    0.2, 0.6, 0.2,
    0.3, 0.3, 0.4
  ), nrow = 3, byrow = TRUE)
  P1_t2 <- matrix(c(
    0.7, 0.2, 0.1,
    0.1, 0.7, 0.2,
    0.2, 0.2, 0.6
  ), nrow = 3, byrow = TRUE)

  P2_t1 <- matrix(c(
    0.6, 0.3, 0.1,
    0.4, 0.5, 0.1,
    0.1, 0.2, 0.7
  ), nrow = 3, byrow = TRUE)
  P2_t2 <- matrix(c(
    0.5, 0.3, 0.2,
    0.2, 0.6, 0.2,
    0.3, 0.3, 0.4
  ), nrow = 3, byrow = TRUE)

  Ps <- list(
    list(P1_t1, P1_t2),
    list(P2_t1, P2_t2)
  )

  out <- summarize_Ps(Ps)

  # Row count: J*(T-1)*K*K = 2*2*3*3 = 36
  expect_equal(nrow(out), 36L)

  # Specific spot checks
  # j = 2, t = 1, from = "1", to = "2" (i.e., row 1, col 2 in P2_t1) = 0.3
  expect_equal(
    out$prob[out$j == 2 & out$t == 1 & out$from == "0" & out$to == "1"],
    0.3
  )
  # j = 1, t = 2, from = "2", to = "3" (row 3, col 3) = 0.6
  expect_equal(
    out$prob[out$j == 1 & out$t == 2 & out$from == "2" & out$to == "2"],
    0.6
  )
})

test_that("summarize_Ps: works for K = 1 (degenerate single-state case)", {
  Ps <- list(list(matrix(1, nrow = 1)))
  out <- summarize_Ps(Ps)

  expect_equal(nrow(out), 1L)
  expect_equal(out$j, 1L)
  expect_equal(out$t, 1L)
  expect_equal(out$from, "0")
  expect_equal(out$to, "0")
  expect_equal(out$prob, 1)
})

test_that("summarize_Ps: no factors are produced for from/to", {
  # Guard against environments where stringsAsFactors might be TRUE
  Ps <- list(list(matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)))
  out <- summarize_Ps(Ps, y_names = c("a", "b"))

  expect_false(is.factor(out$from))
  expect_false(is.factor(out$to))
})

test_that("summarize_Ps: row sums of each matrix are preserved in long form", {
  # Not a strict invariant of the output, but we can check
  # that for each (j, t, from), the to-probabilities sum to the
  # corresponding row sum of the original P matrix.
  P <- matrix(c(0.6, 0.4,
                0.3, 0.7), nrow = 2, byrow = TRUE)
  Ps <- list(list(P))
  out <- summarize_Ps(Ps)

  sum_from0 <- sum(out$prob[out$j == 1 & out$t == 1 & out$from == "0"])
  sum_from1 <- sum(out$prob[out$j == 1 & out$t == 1 & out$from == "1"])

  expect_equal(sum_from0, sum(P[1, ]))
  expect_equal(sum_from1, sum(P[2, ]))
})
