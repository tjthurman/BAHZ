
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
#' A quick little function to turn NaNs and negative numbers into 0
#'
#' DETAILS TO BE ADDED. LINK TO THE GENERAL CLINE EQUATION FUNCTION, MAYBE MOVE
#' THE CLINE PARAMETER DEFINITIONS THERE.
#' @keywords internal
#'
#' @param .df A data frame, containing a column named Fis.
#' Meant to be run as part of a dplyr pipeline
#'
#' @return TO ADD
#'
#' @examples
#' # to be added
#'
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
