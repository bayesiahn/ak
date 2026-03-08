library(ggplot2)
library(dplyr)
library(ggpubr)

y_support <- c(0, 1)

y0_0_dist_initial <- matrix(c(0.75, 0.25), ncol = 1)
y1_1_dist_initial <- matrix(c(0.5, 0.5), ncol = 1)


get_future_dist <- function(f_t, M_t) {
  t(t(f_t) %*% M_t)
}

get_expectation <- function(f_t, y_support) {
  sum(f_t * y_support)
}


get_dist_all <- function(dist_initial, M_first_period) {
  dist_t1 <- get_future_dist(dist_initial, M_first_period)
  return(list(dist_initial, dist_t1))
}

get_plot_dynamics_two_periods_df <- function(y0_0_dist_initial, y1_1_dist_initial,
                   M0_0, M0_1,
                   display_counterfactual_untreated_for_treated = TRUE) {
  # Get the distributions of y for each period
  y0_0_dist_transitions <- get_dist_all(y0_0_dist_initial,M0_0)
  y0_1_dist_transitions <- get_dist_all(y1_1_dist_initial,M0_1)

  y0_0 <- sapply(y0_0_dist_transitions, get_expectation, y_support)
  y0_1 <- sapply(y0_1_dist_transitions, get_expectation, y_support)

  # Create a data frame with x and y coordinates for the lines
  x = rep(1:2, 1)
  y = c(y0_0)
  if (display_counterfactual_untreated_for_treated) {
    x = rep(1:2, 2)
    y = c(y0_1,y0_0)
  }
  df <- data.frame(
    x = x, y = y,
    Outcome = rep(c("Unobserved (Counterfactual)","Observed"),
                  each = 2),
    Treatment_Status = rep(c("Treated","Untreated"),
                           each = 2)
  )

  df
}

plot_dynamics_two_periods <- function(y0_0_dist_initial, y1_1_dist_initial,
                                      M0_0, M0_1,
                                      plot_name = "plots/plot-parallel-trends-violation.pdf",
                                      plot_title = "",
                                      plot_width = 6, plot_height = 6,
                                      display_counterfactual_untreated_for_treated = TRUE,
                                      df_treated = NULL,
                                      draw_arrow_before = FALSE, draw_arrow_after = FALSE) {

  df <- get_plot_dynamics_two_periods_df(y0_0_dist_initial, y1_1_dist_initial,
                                         M0_0, M0_1,
                                         display_counterfactual_untreated_for_treated)
  y_name <- "Average Untreated Potential Outcome"

  if (!is.null(df_treated)) {
    df <- rbind(df, df_treated)
    y_name <- "Average Potential Outcome"
  }


  # Generate a plot
  plot <-
    ggplot(df, aes(x, y, group = interaction(Outcome, Treatment_Status),
                   color = Treatment_Status, linetype = Outcome)) +
    geom_point(aes(x, y, color = Treatment_Status, shape = Treatment_Status),  size = 3) +
    geom_line(aes(group = interaction(Outcome, Treatment_Status))) +
    scale_x_continuous( breaks=c(1,2), limits = c(0.9,2.1)) +
    scale_y_continuous(breaks=seq(0.25,1,0.25), limits = c(0.25, 1)) +
    scale_linetype_manual(values=c("Unobserved (Counterfactual)" = "dashed",
                                   "Observed" = "solid"))+
    scale_color_manual(values=c( "red","blue"))+
    labs(x = "Period", y = y_name, title = "",
         color = "Treatment Cohort",
         shape = "Treatment Cohort",
         linetype = "Outcome")  +
    ggplot2::theme_bw() +
    theme(legend.position='bottom') +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ak::reference_hline(yintercept = 1)

  if (draw_arrow_after) {
    plot <- plot +
      geom_segment(
        x = 2, y = 0.75, xend = 2, yend = 1,
        lineend = "round", linejoin = "round",
        arrow = arrow(length = unit(0.3, "inches")),
        show.legend = FALSE
      )
  }
  if (draw_arrow_before) {
    plot <- plot +  geom_segment(
        x = 1, y = 0.25, xend = 1, yend = 0.5,
        lineend = "round", linejoin = "round",
        arrow = arrow(length = unit(0.3, "inches")),
        show.legend = FALSE
      )
  }

  ggsave(plot_name, plot, width = plot_width, height = plot_height)

  list(df = df, plot = plot)
}

df_treated <- data.frame(
  x = c(1,2),
  y = c(0.5, 0.875),
  Outcome = "Observed",
  Treatment_Status = "Treated"
)


M0_0 <- matrix(c(1/3, 0, 2/3, 1), ncol = 2)
M0_1 <- matrix(c(0, 0, 1, 1), ncol = 2)
slide_plot_width <- 10
slide_plot_height <- 8

plot_parallel_trends <- plot_dynamics_two_periods(y0_0_dist_initial, y1_1_dist_initial, M0_0, M0_1,
  plot_name = "plots/plot-transformation-invariance-parallel-trends.pdf",
  plot_title = "(a) Mean Independence (Parallel Trends)")


plot_parallel_trends_with_treated <- plot_dynamics_two_periods(y0_0_dist_initial, y1_1_dist_initial, M0_0, M0_1,
  plot_name = "plots/plot-transformation-invariance-parallel-trends-with-treated.pdf",
  plot_title = "(a) Mean Independence (Parallel Trends)",
  df_treated = df_treated)


plot_parallel_trends_slides <- plot_dynamics_two_periods(y0_0_dist_initial, y1_1_dist_initial, M0_0, M0_1,
                          plot_name = "plots/plot-slides-parallel-trends-without-counterfactual.pdf",
                          plot_title = "",
                          df_treated = df_treated,
                          plot_width = slide_plot_width, plot_height = slide_plot_height,
                          display_counterfactual_untreated_for_treated = FALSE)
plot_dynamics_two_periods(y0_0_dist_initial, y1_1_dist_initial, M0_0, M0_1,
                          plot_name = "plots/plot-slides-parallel-trends-without-counterfactual-with-arrow.pdf",
                          plot_title = "",
                          df_treated = df_treated,
                          plot_width = slide_plot_width, plot_height = slide_plot_height,
                          display_counterfactual_untreated_for_treated = FALSE,
                          draw_arrow_before = TRUE, draw_arrow_after = FALSE)
plot_dynamics_two_periods(y0_0_dist_initial, y1_1_dist_initial, M0_0, M0_1,
                          plot_name = "plots/plot-slides-parallel-trends-without-counterfactual-with-arrows.pdf",
                          plot_title = "",
                          df_treated = df_treated,
                          plot_width = slide_plot_width, plot_height = slide_plot_height,
                          display_counterfactual_untreated_for_treated = FALSE,
                          draw_arrow_before = TRUE, draw_arrow_after = TRUE)
plot_dynamics_two_periods(y0_0_dist_initial, y1_1_dist_initial, M0_0, M0_1,
                          plot_name = "plots/plot-slides-parallel-trends-with-counterfactual-with-arrows.pdf",
                          plot_title = "",
                          df_treated = df_treated,
                          plot_width = slide_plot_width, plot_height = slide_plot_height,
                          display_counterfactual_untreated_for_treated = TRUE,
                          draw_arrow_before = TRUE, draw_arrow_after = TRUE)
plot_dynamics_two_periods(y0_0_dist_initial, y1_1_dist_initial, M0_0, M0_1,
                          plot_name = "plots/plot-slides-parallel-trends-with-counterfactual.pdf",
                          plot_title = "",
                          df_treated = df_treated,
                          plot_width = slide_plot_width, plot_height = slide_plot_height,
                          display_counterfactual_untreated_for_treated = TRUE,
                          draw_arrow_before = FALSE, draw_arrow_after = FALSE)


M0_0 <- matrix(c(1/3, 0, 2/3, 1), ncol = 2)
M0_1 <- M0_0 # parallel transitions hold

plot_parallel_transitions <- plot_dynamics_two_periods(y0_0_dist_initial, y1_1_dist_initial, M0_0, M0_1,
  plot_name = "plots/plot-transformation-invariance-parallel-transitions.pdf",
  plot_title = "(b) Transition Independence")

plot_parallel_transitions_with_treated <- plot_dynamics_two_periods(y0_0_dist_initial, y1_1_dist_initial, M0_0, M0_1,
  plot_name = "plots/plot-transformation-invariance-parallel-transitions-with-treated.pdf",
  plot_title = "(b) Transition Independence",
  df_treated = df_treated)



plots <- list(plot_parallel_trends$plot, plot_parallel_transitions$plot + labs(y=""))
plots_with_treated <- list(plot_parallel_trends_with_treated$plot,
  plot_parallel_transitions_with_treated$plot + labs(y=""))

plot_combined <- ggpubr::ggarrange(plotlist = plots,
                                  ncol = 2, nrow = 1,
                                  common.legend = TRUE, legend = "bottom")
plot_combined
ggsave("plots/plot-transformation-invariance-combined.pdf", plot_combined,
       width = 8.5, height = 7)


plot_combined_with_treated <- ggpubr::ggarrange(plotlist = plots_with_treated,
                                  ncol = 2, nrow = 1,
                                  common.legend = TRUE, legend = "bottom")
plot_combined_with_treated
ggsave("plots/plot-transformation-invariance-combined-with-treated.pdf", plot_combined_with_treated,
       width = 8.5, height = 7)

