
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

#' Parse the prior yaml file and check it
#'
#' DESCRIPTION TO BE ADDED
#'
#' Used interanlly, in \link{create_cline_model} and \link{make_init_list}.
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



# Extract values out of priors ----------------------------------------------

#
#' Extract the first value from comma-separated list/vector
#'
#' Uses regular expressions implmented in stringr package to extract the first
#' value, that is, the one after an open parenthesis and before a comma.
#' Properly handles whitespace and decimals.
#'
#' Used interanlly, in \link{init_single_chain}
#'
#'
#' @keywords internal
#'
#' @param string The string to extract a value from
#'
#' @return A numeric value
#'
#'

extract_first <- function(string) {
  res <- stringr::str_extract(string, "\\([:blank:]*[0-9]*\\.*[0-9]*[:blank:]*,") %>%
    stringr::str_remove_all("[(,]") %>%
    stringr::str_squish() %>%
    as.numeric
  assertthat::assert_that(is.na(res) == F, msg = "Could not parse prior. Check your prior file!")
  res
}

#
#' Extract the last value from comma-separated list/vector
#'
#' Uses regular expressions implmented in stringr package to extract the last
#' value, that is, the one after a comma and before a close parenthesis.
#' Properly handles whitespace and decimals.
#'
#' Used interanlly, in \link{init_single_chain}
#'
#'
#' @keywords internal
#'
#' @param string The string to extract a value from
#'
#' @return A numeric value
#'
#'

extract_last <- function(string) {
  res <- stringr::str_extract(string, ",[:blank:]*[0-9]*\\.*[0-9]*[:blank:]*\\)") %>%
    stringr::str_remove_all("[,)]") %>%
    stringr::str_squish() %>%
    as.numeric
  assertthat::assert_that(is.na(res) == F, msg = "Could not parse prior. Check your prior file!")
  res
}

#
#' Extract the value from a single-parameter distribution
#'
#' Uses regular expressions implmented in stringr package to extract the only
#' value from a distribution with one parameter. That is, is extracts the
#' numbers between two parentheses. Properly handles whitespace and decimals.
#'
#' Used interanlly, in \link{init_single_chain}
#'
#'
#' @keywords internal
#'
#' @param string The string to extract a value from
#'
#' @return A numeric value
#'
#'
extract_only <- function(string) {
  res <- stringr::str_extract(string, "\\([:blank:]*[0-9]*\\.*[0-9]*[:blank:]*\\)") %>%
    stringr::str_remove_all("[,)(]") %>%
    stringr::str_squish() %>%
    as.numeric
  assertthat::assert_that(is.na(res) == F, msg = "Could not parse prior. Check your prior file!")
  res
}

