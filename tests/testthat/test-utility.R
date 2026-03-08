test_that("y_to_y_indices returns correct results for vectors", {
  # Setup
  y <- c("a", "b", "c")
  y_values <- c("a", "b", "c", "d")
  expected_result <- c(1, 2, 3)

  # Exercise
  result <- y_to_y_indices(y, y_values)

  # Verify
  expect_equal(as.numeric(result), expected_result)
})

test_that("y_to_y_indices returns correct results for matrices", {
  # Setup
  y <- matrix(c("a", "b", "c", "d", "e", "f"), nrow = 2)
  y_values <- c("a", "b", "c", "d", "e", "f")
  expected_result <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)

  # Exercise
  result <- y_to_y_indices(y, y_values)

  # Verify
  expect_equal(result, expected_result)
})

test_that("y_to_y_indices handles edge cases", {
  # Setup
  y <- c("a", "b", "c")
  y_values <- c()
  expected_result <- NA

  # Exercise
  result <- y_to_y_indices(y, y_values)

  # Verify
  expect_equal(result, expected_result)
})

test_that("transitions_to_transition_indicators returns correct dimensions", {
  # [,1] [,2] [,3]
  # [1,]    1    2    2
  # [2,]    2    1    2
  y_indices_matrix <- matrix(c(1, 2, 2, 1, 2, 2), nrow = 2)
  y_values <- 1:2
  result <- transitions_to_transition_indicators(y_indices_matrix, y_values)

  # Check dimensions
  expect_equal(length(result), ncol(y_indices_matrix) - 1)
  expect_equal(length(result[[1]]), nrow(y_indices_matrix))
  expect_equal(dim(result[[1]][[1]]), c(length(y_values), length(y_values)))

  # Check consistency
  expect_equal(result[[1]][[1]][1,2], 1)
  expect_equal(result[[2]][[1]][2,2], 1)
  expect_equal(result[[1]][[2]][2,1], 1)
  expect_equal(result[[2]][[2]][1,2], 1)
})

test_that("transitions_to_transition_indicators returns correct transition matrices", {
  # Case 1
  # [,1] [,2] [,3]
  # [1,]    1    2    2
  # [2,]    2    1    2
  y_indices_matrix <- matrix(c(1, 2, 2, 1, 2, 2), nrow = 2)
  y_values <- 1:2
  result <- transitions_to_transition_indicators(y_indices_matrix, y_values)

  expect_equal(result[[1]][[1]][1,2], 1)
  expect_equal(result[[2]][[1]][2,2], 1)
  expect_equal(result[[1]][[2]][2,1], 1)
  expect_equal(result[[2]][[2]][1,2], 1)

  # Fixing current period as initial period
  result <- transitions_to_transition_indicators(y_indices_matrix, y_values, 1)

  expect_equal(result[[1]][[1]][1,2], 1)
  expect_equal(result[[2]][[1]][1,2], 1)
  expect_equal(result[[1]][[2]][2,1], 1)
  expect_equal(result[[2]][[2]][2,2], 1)

  # Case 2
  # [,1] [,2]
  # [1,]    1    2
  # [2,]    2    1
  y_indices_matrix <- matrix(c(1, 2, 2, 1), nrow = 2)
  y_values <- 1:2
  expected_result <- list(
    list(matrix(c(0, 0, 1, 0), nrow = 2), matrix(c(0, 1, 0, 0), nrow = 2))
  )
  result <- transitions_to_transition_indicators(y_indices_matrix, y_values)
  result_0 <- transitions_to_transition_indicators(y_indices_matrix, y_values, 1)
  expect_equal(result, expected_result)
  expect_equal(result_0, expected_result)
})

test_that("transitions_to_transition_indicators handles edge cases", {
  # Edge Case: Single row matrix
  y_indices_matrix <- matrix(1:3, nrow = 1)
  y_values <- 1:3
  result <- transitions_to_transition_indicators(y_indices_matrix, y_values)
  expect_equal(length(result), length(y_indices_matrix) - 1)
  expect_equal(dim(result[[1]][[1]]), c(length(y_values), length(y_values)))

  # Edge Case: Single column matrix
  y_indices_matrix <- matrix(1:2, ncol = 1)
  y_values <- 1:2
  result <- transitions_to_transition_indicators(y_indices_matrix, y_values)
  expect_equal(result, list())
})

test_that("y_indices_to_initial_outcome_indicators returns correct one-hot encoded matrix", {
  y_indices_matrix <- matrix(c(1, 2, 3, 1, 3, 2), nrow = 3)
  y_values <- 1:3
  expected_result <- matrix(c(1,0,0, 0,1,0, 0,0,1), nrow = 3, byrow = TRUE)
  result <- y_indices_to_initial_outcome_indicators(y_indices_matrix, y_values)

  expect_equal(result, expected_result)
})

test_that("y_indices_to_initial_outcome_indicators handles edge cases", {
  # Edge Case: Empty matrix
  y_indices_matrix <- matrix(numeric(0), nrow = 0)
  y_values <- numeric(0)
  result <- y_indices_to_initial_outcome_indicators(y_indices_matrix, y_values)
  expect_equal(dim(result), c(0, length(y_values)))

  # Edge Case: Single row matrix
  y_indices_matrix <- matrix(1:3, nrow = 1)
  y_values <- 1:3
  result <- y_indices_to_initial_outcome_indicators(y_indices_matrix, y_values)
  expect_equal(dim(result), c(1, length(y_values)))

  # Edge Case: Single column matrix (only initial states)
  y_indices_matrix <- matrix(1:2, ncol = 1)
  y_values <- 1:2
  result <- y_indices_to_initial_outcome_indicators(y_indices_matrix, y_values)
  expect_equal(dim(result), c(2, length(y_values)))
})

test_that("y_indices_to_initial_outcome_indicators ensures row sums are one", {
  y_indices_matrix <- matrix(c(1, 2, 3, 1, 3, 2), nrow = 3)
  y_values <- 1:3
  result <- y_indices_to_initial_outcome_indicators(y_indices_matrix, y_values)

  expect_true(all(rowSums(result) == 1))
})


test_that("multiply_matrices_in_list returns correct product", {
  K_max <- 3

  for (K in 1:K_max) {
    P_1 <- matrix(runif(K^2), nrow = K)
    P_2 <- matrix(runif(K^2), nrow = K)
    P_3 <- matrix(runif(K^2), nrow = K)


    matrices <- list(P_1)
    matrix <- multiply_matrices_in_list(matrices)
    expect_true(all(matrix == P_1))

    matrices <- list(P_1, P_2)
    matrix <- multiply_matrices_in_list(matrices)
    expect_true(all(matrix == P_1 %*% P_2))

    matrices <- list(P_1, P_2, P_3)
    matrix <- multiply_matrices_in_list(matrices)
    expect_true(all(matrix == P_1 %*% P_2 %*% P_3))
  }
})

test_that("concatenate_Ps correctly concatenates matrices", {
  Ps_pre <- list(list(matrix(1:4, nrow = 2)), list(matrix(5:8, nrow = 2)))
  Ps_post <- list(list(matrix(9:12, nrow = 2)), list(matrix(13:16, nrow = 2)))
  expected_Ps <- list(list(matrix(1:4, nrow = 2), matrix(9:12, nrow = 2)),
                      list(matrix(5:8, nrow = 2), matrix(13:16, nrow = 2)))

  result <- concatenate_Ps(Ps_pre, Ps_post)

  expect_equal(length(result), length(expected_Ps))
  for (j in 1:length(result)) {
    expect_equal(result[[j]], expected_Ps[[j]])
  }
})

test_that("concatenate_Ps handles mismatched list lengths", {
  Ps_pre <- list(matrix(1:4, nrow = 2))
  Ps_post <- list(matrix(9:12, nrow = 2), matrix(13:16, nrow = 2))

  expect_error(concatenate_Ps(Ps_pre, Ps_post))
})

test_that("long_df_to_wide returns correct wide format matrix", {
  data <- data.frame(id = c(1,1,1,2,2,2), t = c(1,2,3,1,2,3),
                     y = c(1,2,3,4,5,6), g = c(1,1,1,2,2,2))
  result <- long_df_to_wide(data, "y", "id", "t")
  expected <- matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE)
  expect_equal(result, expected, ignore_attr = TRUE)
})

test_that("long_df_to_y_and_g extracts y matrix and g vector correctly", {
  data <- data.frame(id = c(1,1,1,2,2,2), t = c(1,2,3,1,2,3),
                     y = c(1,2,3,4,5,6), g = c(1,1,1,0,0,0))
  result <- long_df_to_y_and_g(data, "y", "g", "id", "t")
  expected_y <- matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE)
  expected_g <- c(1,0)
  expect_equal(result$y, expected_y, ignore_attr = TRUE)
  expect_equal(result$g, expected_g, ignore_attr = TRUE)
})

test_that("long_df_to_y_and_g handles missing time points correctly", {
  data <- data.frame(id = c(1,1,2,2), t = c(1,3,1,3),
                     y = c(1,3,4,6), g = c(1,1,3,3))
  result <- long_df_to_y_and_g(data, "y", "g", "id", "t")
  expected_y <- matrix(c(1,3,4,6), nrow = 2, byrow = TRUE)
  expected_g <- c(1,2) # c(1,2) because 3 is rescaled to 2
  expect_equal(result$y, expected_y, ignore_attr = TRUE)
  expect_equal(result$g, expected_g, ignore_attr = TRUE)
})

