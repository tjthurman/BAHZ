
# Utility functions for cap -----------------------------------------------
# Not user facing, just used internally.


# Pipe operator -----------------------------------------------------------
#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


# Correct_Fis -------------------------------------------------------------
#
#' Correct Fis estimates
#'
#' A function that corrects Fis estimates which are negative or NaN, turns them
#' to 0. Reads in a data frame with an Fis column. Meant to be part of a dplyr
#' pipeline.
#'
#' Used interanlly, in \link{sim_data_from_cline}.
#'
#' @keywords internal
#'
#' @param .df A data frame, containing a column named Fis (this is not checked).
#'
#' @return The supplied data frame, with a corrected Fis column
#'
#'

correct_fis <- function(.df) {
  for (element in 1:length(.df$Fis)) {
    if (is.nan(.df$Fis[element]) == T) {
      .df$Fis[element] <- 0
    }
    else if (.df$Fis[element] < 0) {
      .df$Fis[element] <- 0
    }
  }
 .df
}
