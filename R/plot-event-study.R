#' Plot an Event-Study from a DID Result
#'
#' Produces a clean event-study plot from a `did` estimation object, using the
#' aggregated dynamic treatment effects obtained via `did::aggte()`.
#' The function adjusts the post-event period labels to start from 1 (rather than 0)
#' and formats the plot with a minimalist publication-style theme.
#'
#' @param did_result An object returned by `did::att_gt()`, representing
#'   group-time average treatment effect estimates.
#' @param uniform_ci Logical; if `TRUE`, plots uniform (simultaneous) confidence
#'   intervals instead of pointwise ones. Passed to the `cband` argument in
#'   `did::aggte()`. Defaults to `FALSE`.
#'
#' @return A `ggplot` object showing event-time dynamic ATT estimates with
#'   confidence intervals, using:
#'   - integer x-axis labels (post-event periods start at 1),
#'   - no minor grid lines on the x-axis,
#'   - customized axis labels and no title,
#'   - legend placed at the bottom.
#'
#' @details
#' The function calls:
#' \enumerate{
#'   \item `did::aggte()` with \code{type = "dynamic"} to compute event-study estimates;
#'   \item `did::ggdid()` to generate the base event-study plot;
#'   \item adjusts x-axis labels so that the post-event periods are labeled starting from 1;
#'   \item applies consistent theming for clean presentation.
#' }
#'
#' @examples
#' \dontrun{
#' # Example using the simulated dataset from the did package
#' data(simdata)
#' attgt <- did::att_gt(yname = "Y",
#'                      tname = "year",
#'                      idname = "id",
#'                      gname = "G",
#'                      xformla = ~ X,
#'                      data = simdata)
#' plot_did_result(attgt, uniform_ci = TRUE)
#' }
#'
#' @export
plot_did_result <- function(did_result, uniform_ci = FALSE) {
  p <- did::ggdid(did::aggte(did_result,
                             type = "dynamic",
                             cband = uniform_ci))

  p$layers <- p$layers[!sapply(p$layers, function(l) {
    inherits(l$geom, "GeomVline") || inherits(l$geom, "GeomHline")
  })]

  suppressMessages(
    p <- p +
    reference_hline() +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "dashed") +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(c(0, ggplot_build(p)$data[[1]]$x))),
      labels = function(x) x + 1
    ) +
    ggplot2::scale_color_manual(values = rep("black", 2)) +
    ggplot2::scale_fill_manual(values = rep("black", 2)) +
    ggplot2::guides(color = "none", fill = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::theme(panel.grid.minor.x = element_blank()) +
    ggplot2::ylab(did_result$empirical_spec$y_lab_average_1) +
    ggplot2::xlab(did_result$empirical_spec$x_lab) +
    ggplot2::ggtitle(NULL)
  )
  p
}
