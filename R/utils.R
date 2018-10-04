
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

# Parse prior file --------------------------------------------------------
#
#' Parse the prior yaml file and check it
#'
#' DESCRIPTION TO BE ADDED
#'
#' Used interanlly, in LINK TO FUNCTIONS THAT USE IT. REMOVE EXPORT WHEN DONE WITH TESTING.
#'
#' @export
#'
#' @keywords internal
#'
#' @param prior_file filepath to the prior file
#'
#' @return a named list containing the specifed priors
#'
#'

parse_prior_file <- function(prior_file){
  priors <- yaml::yaml.load_file(prior_file, as.named.list = T)
  # ADD a check to make sure it is 11 long, and all the right names are there

  assertthat::assert_that(length(priors) == 11, msg = "Incorrect number of priors, there should be 11!\nDouble-check your prior file")
  name.check <- names(priors) %in% c("center", "width", "pmin", "pmax",
                                     "deltaL", "deltaR", "deltaM",
                                     "tauL", "tauR", "tauM", "f")

  if (sum(name.check) != 11) {
    offenders <- as.vector(names(priors)[which(name.check == F)])
    stop(paste("\n", toString(offenders), "\nis/are not valid parameter names.\nDouble-check your prior file", sep = ""))
  }
  priors
}
