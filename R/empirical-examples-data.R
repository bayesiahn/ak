#' @title Data Transformation Functions for Empirical Examples
#'
#' @description
#' Functions for transforming model outputs into data frames suitable for
#' analysis and visualization in empirical example analyses.
#'
#' @name empirical-examples-data
NULL

#' Generate Long-Format Data Frame from Model
#'
#' Extracts outcome and group data from a model object and transforms it into
#' a long-format data frame suitable for analysis and visualization.
#'
#' @param model A model object containing outcome and group information.
#' @return A data frame in long format with columns for id, time period (t),
#'   outcome value (y), and group identifier (g).
#' @export
generate_df_y <- function(model) {
  y_g <- get_y_numeric_and_g_from_model(model)
  y <- y_g$y
  g <- y_g$g

  df_y <- y %>% as.data.frame()
  rownames(df_y) <- paste0(1:nrow(df_y))
  colnames(df_y) <- paste0(1:ncol(df_y))

  df_g <- data.frame(id = as.character(1:length(g)), g = g)
  df_y %>%
    tibble::rownames_to_column("id") %>%
    tidyr::pivot_longer(cols = -id,
                        names_to = "t",
                        values_to = "y") %>%
    dplyr::left_join(df_g, by = "id")
}

#' Generate Pretreatment Outcome Data
#'
#' Extracts outcome values at the pretreatment period (one period before
#' treatment) for each observation in the model.
#'
#' @param model A model object containing outcome data and treatment period
#'   information.
#' @return A data frame with columns for id and y_at_pretreatment.
#' @export
generate_df_y_pretreatment <- function(model) {
  treated_period <- model$treated_period

  df_y_pretreatment <- generate_df_y(model) %>%
    dplyr::filter(t == treated_period - 1) %>%
    dplyr::mutate(y_at_pretreatment = y) %>%
    dplyr::select(id, y_at_pretreatment)
}
