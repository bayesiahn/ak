test_that("model_to_scores returns correct structure", {
  set.seed(123)
  N <- 100
  T_max <- 4
  J <- 2
  g <- sample(c(0, 1), N, replace = TRUE)
  y <- matrix(sample(1:3, N * T_max, replace = TRUE), nrow = N, ncol = T_max)
  model <- generate_transition_model(outcomes = 1:3, J = J, T_max = T_max)

  scores <- model_to_scores(y, g, model)

  expect_type(scores, "list")
  expect_true(all(c("Y", "treated_period", "priors", "priors_treated",
                    "Ps_control", "Ps_treated", "Ps_control_empirical",
                    "Ps_treated_empirical", "pmfs_initial_control",
                    "pmfs_initial_treated", "pmfs_pre_treatment_control",
                    "pmfs_pre_treatment_treated") %in% names(scores)))

  expect_failure(expect_true("p_y1d" %in% names(scores)))
})

# test_that("scores_to_model returns correct class and structure", {
#   set.seed(123)
#   model <- generate_transition_model(outcomes = 1:2, J = 2, T_max = 2)
#   sample <- generate_sample(model, N = 50000)
#   y <- sample$y
#   g <- sample$g
#   scores <- model_to_scores(y, g, model)
#   weights <- sample(c(-1, 1), nrow(y), replace = TRUE)
#   model_from_scores <- scores_to_model(scores, weights)
#
#   # Check if multiplier bootstrap estimates are centered around the actual params
#   atts <- vector()
#   ltatts <- list()
#   for (k in 1:100) {
#     weights <- sample(c(-1, 1), nrow(y), replace = TRUE)
#     model_from_scores <- scores_to_model(scores, weights)
#     summary <- summarize_transition_model(model_from_scores)
#     atts <- c(atts, summary$ATT)
#     ltatts[[k]] <- unlist(summary$LTATTs)
#   }
#   atts_mean <- median(atts)
#   ltatts_mean <- apply(do.call(rbind, ltatts), 2, median)
#
#   summary <- summarize_transition_model(model)
#   att_actual <- summary$ATT
#   ltatts_actual <- unlist(summary$LTATTs)
#   expect_equal(atts_mean, att_actual, tolerance = 1e0)
#   expect_equal(ltatts_mean, ltatts_actual, tolerance = 1e0)
#
#   # Check if contains all elements
#   expect_type(model_from_scores, "list")
#   expect_equal(class(model_from_scores), "TransitionModel")
#   expect_true(all(c("priors", "Ps_control", "pmfs_initial_control") %in% names(model_from_scores)))
# })
