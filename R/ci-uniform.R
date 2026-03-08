#' @title Confidence intervals from bootstrap paths (uniform or pointwise)
#' @description Construct confidence intervals for a sequence of estimates
#'   observed over periods \eqn{t = 1, \ldots, T} using bootstrap studentization.
#'   When \code{uniform_ci = TRUE}, builds uniform (simultaneous) bands using
#'   the supremum of t-statistics across periods. When \code{uniform_ci = FALSE},
#'   builds pointwise bands using per-period critical values.
#' @param estimate A numeric vector of length \eqn{T} containing the point
#'   estimates. If unnamed, indices \code{"1"}, \code{"2"}, etc. are used as
#'   period labels.
#' @param bootstrap_estimates A list of length \eqn{B}; each element is a
#'   numeric vector of length \eqn{T} with bootstrap estimates for all periods.
#' @param uniform_ci A logical; if \code{TRUE}, compute uniform (simultaneous)
#'   confidence intervals; if \code{FALSE}, compute pointwise intervals.
#' @return A \code{data.frame} with columns \code{t} (period labels),
#'   \code{estimate} (point estimate), \code{sd} (bootstrap standard deviation),
#'   \code{critical_value}, \code{ci_lower}, and \code{ci_upper}.
#' @examples
#' est <- c(a = 1.0, b = 1.2, c = 0.9)
#' boot <- list(
#'   c(0.9, 1.1, 1.0),
#'   c(1.1, 1.3, 0.8),
#'   c(1.0, 1.2, 0.9)
#' )
#' # Pointwise (default)
#' bootstrap_estimates_to_ci_df(est, boot)
#' # Uniform (simultaneous)
#' bootstrap_estimates_to_ci_df(est, boot, uniform_ci = TRUE)
#' @importFrom stats sd quantile
#' @export
bootstrap_estimates_to_ci_df <- function(estimate, bootstrap_estimates,
                                         uniform_ci = FALSE) {
  # ---- Input checks -------------------------------------------------------
  if (!is.numeric(estimate) || length(estimate) == 0L) {
    stop("`estimate` must be a non-empty numeric vector.")
  }
  if (!is.list(bootstrap_estimates) || length(bootstrap_estimates) == 0L) {
    stop("`bootstrap_estimates` must be a non-empty list of numeric vectors.")
  }
  T_len <- length(estimate)

  # Ensure every bootstrap vector matches length T
  lens <- vapply(bootstrap_estimates, length, integer(1))
  if (any(lens != T_len)) {
    stop("All elements of `bootstrap_estimates` must have the same length as `estimate`.")
  }

  # Period labels t: use names if available, else "1","2",...,"T"
  t_lab <- names(estimate)
  if (is.null(t_lab) || is.na(t_lab[1])) {
    t_lab <- as.character(seq_len(T_len))
    names(estimate) <- t_lab
  } else {
    t_lab <- as.character(t_lab)
  }

  # ---- Arrange bootstrap into B x T matrix --------------------------------
  B <- length(bootstrap_estimates)
  boot_mat <- do.call(rbind, lapply(bootstrap_estimates, as.numeric))
  boot_mat <- as.matrix(boot_mat)
  colnames(boot_mat) <- t_lab

  # ---- 1) Bootstrap SD by period ------------------------------------------
  sd_vec <- apply(boot_mat, 2L, stats::sd)

  # ---- 2) Bootstrap t-stats for each (b,t) --------------------------------
  est_vec <- as.numeric(estimate)
  tstat_mat <- sweep(boot_mat, 2L, est_vec, FUN = "-")
  zero_sd <- sd_vec == 0 | is.na(sd_vec)
  if (any(!zero_sd)) {
    tstat_mat[, !zero_sd] <- sweep(tstat_mat[, !zero_sd, drop = FALSE],
                                   2L, sd_vec[!zero_sd], FUN = "/")
  }
  if (any(zero_sd)) {
    tstat_mat[, zero_sd] <- 0
  }

  # ---- 3) Critical value(s) -----------------------------------------------
  if (isTRUE(uniform_ci)) {
    # sup |t| per bootstrap and its 0.95 quantile (scalar)
    max_t_by_b <- apply(abs(tstat_mat), 1L, max, na.rm = TRUE)
    crit <- rep(as.numeric(stats::quantile(max_t_by_b, probs = 0.95,
                                           na.rm = TRUE, type = 7)), T_len)
  } else {
    # pointwise: for each t, 0.95 quantile of |t| distribution (vector length T)
    abs_t <- abs(tstat_mat)
    crit <- apply(abs_t, 2L, function(col)
      as.numeric(stats::quantile(col, probs = 0.95, na.rm = TRUE, type = 7)))
  }

  # ---- 4) CIs --------------------------------------------------------------
  ci_half_width <- crit * sd_vec
  ci_lower <- est_vec - ci_half_width
  ci_upper <- est_vec + ci_half_width

  # ---- Output data.frame ---------------------------------------------------
  out <- data.frame(
    t = t_lab,
    estimate = as.numeric(estimate),
    sd = as.numeric(sd_vec),
    critical_value = as.numeric(crit),
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    stringsAsFactors = FALSE
  )
  out
}
