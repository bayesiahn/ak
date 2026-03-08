#' Summarize a list of transition matrices
#'
#' This function summarizes a list of transition matrices by converting them
#' into a more interpretable format.
#' ≠
#' @param Ps A J-list of (T-1)-list of transition matrices to be summarized.
#' Each transition matrix is a KxK matrix representing the probabilities of
#' transitioning from one state (row) to another (column).
#' @param y_names (Optional) A character vector of names for the outcome states.
#' If not provided, default names will be used.
#' @return A data frame object whose columns `j` `t` `from` `to` `prob` represent
#' the transition probabilities for each state in each period for each latent class.
#' @export
summarize_Ps <- function(Ps, y_names = NULL) {
  J <- length(Ps)
  Tm1 <- length(Ps[[1]])
  K <- nrow(Ps[[1]][[1]])
  if (is.null(y_names)) y_names <- as.character(0:(K-1))
  stopifnot(length(y_names) == K)

  # fixed (from,to) grid; we'll reuse it for every (j,t)
  grid <- data.frame(
    from = y_names[rep(seq_len(K), each = K)],
    to   = y_names[rep(seq_len(K), times = K)],
    stringsAsFactors = FALSE
  )

  out <- vector("list", J * Tm1)
  idx <- 1L
  for (j in seq_len(J)) {
    for (t in seq_len(Tm1)) {
      P <- Ps[[j]][[t]]
      # clone grid and add prob/j/t
      df <- grid
      df$prob <- as.vector(t(P))  # row-major to match (from,to)
      df$j <- j
      df$t <- t
      # order columns once here
      out[[idx]] <- df[, c("j", "t", "from", "to", "prob")]
      idx <- idx + 1L
    }
  }
  do.call(rbind, out)
}
