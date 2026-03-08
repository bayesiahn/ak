#' Print method for empirical_analysis_result objects
#'
#' Summarizes an empirical analysis result, including transition model and DiD outputs.
#'
#' @param x An object of class \code{"empirical_analysis_result"}.
#' @param ... Additional arguments (ignored).
#' @export
print.empirical_analysis_result <- function(x, ...) {
  cat("Empirical Analysis Result\n")
  cat("--------------------------\n")

  # 1. Specification Summary
  if (!is.null(x$empirical_spec)) {
    print(x$empirical_spec)
  }

  # 2. Transition Model Summary
  if (!is.null(x$transition_models)) {
    n_tm <- length(x$transition_models)
    cat("\nTransition Models Estimated: ", n_tm, "\n")

    att_vec <- numeric(n_tm)

    for (j in seq_len(n_tm)) {
      tm <- x$transition_models[[j]]
      tm_summary <- summarize_transition_model(tm)
      att_vec[j] <- mean(tm_summary$ATT)

      cat("\nTransition Model (J = ", j, "):\n", sep = "")
      cat("  Aggregate ATT: ", round(mean(tm_summary$ATT), 3), "\n")
      cat("  Model AIC: ", round(mean(tm_summary$aic), 3), "\n")
      cat("  Model BIC: ", round(mean(tm_summary$bic), 3), "\n")

      # Group shares
      priors <- tm$priors
      priors_treated <- tm$priors_treated
      group_names <- paste0("Group ", seq_along(priors))
      cat("  Latent Type Shares (entire population):\n")
      print(setNames(round(priors, 3), group_names))
      cat("  Latent Type Shares (treated population):\n")
      print(setNames(round(priors_treated, 3), group_names))

      # LTATT matrix
      LTATTs_matrix <- tryCatch(
        do.call(cbind, tm_summary$LTATTs),
        error = function(e) NULL
      )
      if (!is.null(LTATTs_matrix)) {
        rownames(LTATTs_matrix) <- paste0("t = ", seq_len(nrow(LTATTs_matrix)))
        colnames(LTATTs_matrix) <- group_names

        cat("  Latent Type ATT (LTATTs by time):\n")
        print(round(LTATTs_matrix, 3))

        # Aggregate LTATT per group
        agg_ltatt <- colMeans(LTATTs_matrix, na.rm = TRUE)
        cat("  Aggregate LTATT per Latent Type:\n")
        print(round(agg_ltatt, 3))
      }
    }

    cat("\nAggregate ATT (by model):\n")
    att_named <- setNames(round(att_vec, 3), paste0("J = ", seq_len(n_tm)))
    print(att_named)
  }

  # 3. DiD Summary
  if (!is.null(x$did_result)) {
    cat("\nDifference-in-Differences (DiD) Estimate:\n")
    treated_period <- max(x$did_result$group)
    T_max <- length(x$did_result$att) + 1

    att_aggregate <- x$did_result$att_aggregate
    att_vec <- x$did_result$att[(treated_period - 1):(T_max - 1)]
    se_vec <- x$did_result$se[(treated_period - 1):(T_max - 1)]
    post_treat_periods <- paste0("t = ", treated_period:T_max)

    cat("  Aggregate ATT: ", formatC(att_aggregate, digits = 3, format = "f"), "\n")
    cat("  ATT Estimates by Time Period (standard errors in parenthesis):\n")
    for (i in seq_along(att_vec)) {
      cat("    ", post_treat_periods[i], ": ",
          formatC(att_vec[i], digits = 3, format = "f"),
          " (", formatC(se_vec[i], digits = 3, format = "f"), ")\n", sep = "")
    }
  }

  invisible(x)
}
