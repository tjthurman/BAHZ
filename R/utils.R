# Utility functions for bahz -----------------------------------------------
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
#' Used internally, in \code{\link{sim_geno_cline}}.
#'
#' @keywords internal
#'
#' @param .df A data frame, containing a column named Fis.
#'
#' @return The supplied data frame, with a corrected Fis column
#'
#'

correct_fis <- function(.df) {
  assertthat::assert_that(("Fis" %in% names(.df)) == T,
                          msg = "Can't correct Fis, no column named Fis in dataframe")
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
#' Reads in the yaml file containing the prior, checking for the proper number
#' and names of priors.
#'
#' Used internally, in \code{\link{prep_prior_list}} and
#' \code{\link{prep_init_list}}.
#'
#'
#' @keywords internal
#'
#'
#' @param prior_file filepath to the prior file
#'
#' @return a named list containing the specifed priors
#'
#'

parse_prior_file <- function(prior_file) {
  path_to_prior <-
    file.path(normalizePath(prior_file), fsep = .Platform$file.sep)
  priors <- yaml::yaml.load_file(path_to_prior, as.named.list = T)
  # ADD a check to make sure it is 11 long, and all the right names are there

  assertthat::assert_that(length(priors) == 11, msg = "Incorrect number of priors, there should be 11!\nDouble-check your prior file")
  name.check <-
    names(priors) == c("center", "width",
                       "pmin",  "pmax",
                       "deltaL", "deltaR",
                       "deltaM", "tauL",
                       "tauR",  "tauM",
                       "f")

  if (sum(name.check) != 11) {
    offenders <- as.vector(names(priors)[which(name.check == F)])
    stop(
      paste(
        "\n",
        toString(offenders),
        "\nis/are either invalid parameter names or out of order.\nDouble-check your prior file against the provided template",
        sep = ""
      )
    )
  }
  priors
}



# Extract values out of priors ----------------------------------------------


#' Extract values for prior distributions
#'
#' @name extractValue
#'
#' @description Internal functions used in \code{\link{init_single_chain}} to extract
#' numerical values from priors. All use regular expressions as implemented
#' in the \code{\link{stringr}} package, and all properly handle whitespace and
#' decimals.
#'
#' \code{extract_first} gets the first value, that is, the one after an open
#' parenthesis and before a comma.
#'
#' \code{extract_last} gets the last value, that is, the value after a comma and
#' before a close parenthesis.
#'
#' \code{extract_only} gets the only value from a distribution with one
#' parameter, that is, it gets the numbers between two parentheses.
#'
#'
#' @keywords internal
#'
#' @param string The string to extract a value from
#'
#' @return A numeric value
#'

NULL


#'
#' Extract values for prior distributions
#'
#' @rdname extractValue
#'
#'

extract_first <- function(string) {
  assertthat::assert_that(is.character(string) == T, msg = "Could not parse prior. Check your prior file!")
  res <-
    stringr::str_extract(string, "\\([:blank:]*[0-9]*\\.*[0-9]*[:blank:]*,") %>%
    stringr::str_remove_all("[(,]") %>%
    stringr::str_squish() %>%
    as.numeric
  assertthat::assert_that(length(res) == 1, msg = "Could not parse prior. X Check your prior file!")
  assertthat::assert_that(is.na(res) == F, msg = "Could not parse prior. Check your prior file!")
  res
}


#' @rdname extractValue
#'
#'

extract_last <- function(string) {
  assertthat::assert_that(is.character(string) == T, msg = "Could not parse prior. Check your prior file!")
  res <-
    stringr::str_extract(string, ",[:blank:]*[0-9]*\\.*[0-9]*[:blank:]*\\)") %>%
    stringr::str_remove_all("[,)]") %>%
    stringr::str_squish() %>%
    as.numeric
  assertthat::assert_that(length(res) == 1, msg = "Could not parse prior. Check your prior file!")
  assertthat::assert_that(is.na(res) == F, msg = "Could not parse prior. Check your prior file!")
  res
}


#'
#' Extract values for prior distributions
#'
#' @rdname extractValue
#'
#'
extract_only <- function(string) {
  assertthat::assert_that(is.character(string) == T, msg = "Could not parse prior. Check your prior file!")
  res <-
    stringr::str_extract(string, "\\([:blank:]*[0-9]*\\.*[0-9]*[:blank:]*\\)") %>%
    stringr::str_remove_all("[,)(]") %>%
    stringr::str_squish() %>%
    as.numeric
  assertthat::assert_that(length(res) == 1, msg = "Could not parse prior. Check your prior file!")
  assertthat::assert_that(is.na(res) == F, msg = "Could not parse prior. Check your prior file!")
  res
}


# Check that priors are supported ---------------------------------------

#' Check that the distribution specified for a given parameter is supported
#'
#'
#' Used internally, in \code{\link{prep_prior_list}}.
#'
#' @keywords internal
#'
#'
#' @param parameter the parameter for which a prior distribution is specified
#'
#' @param distribution the prior distribution for that parameter
#'
#' @return TRUE/FALSE: is the prior supported for that parameter?
#'


check_prior_supported <- function(parameter, distribution) {
  result <- FALSE
  if (parameter %in% c("center", "width")) {
    if (distribution %in% c("normal", "uniform")) {
      result <- T
    }
  }
  if (parameter %in% c("pmin", "pmax")) {
    if (distribution %in% c("uniform")) {
      result <- T
    }
  }
  if (parameter %in% c("deltaL", "deltaR", "deltaM")) {
    if (distribution %in% c("exponential")) {
      result <- T
    }
  }
  if (parameter %in% c("tauL", "tauR", "tauM")) {
    if (distribution %in% c("uniform")) {
      result <- T
    }
  }
  if (parameter %in% c("f")) {
    if (distribution %in% c("uniform")) {
      result <- T
    }
  }
  result
}

# Count commas for checking prior specification --------------------------


#' Check that the distribution specified has the right number of parameters
#'
#'
#' Used internally, in \code{\link{prep_prior_list}}.
#'
#' @keywords internal
#'
#' @param distribution the prior distribution
#'
#' @param string the string specifiying the parameter values, parsed from the prior_config file
#'
#' @return TRUE/FALSE: is the prior specified correctly?
#'
check_prior_specification <- function(distribution, string) {
  result <- F
  num.commas <- stringr::str_count(string, "\\,")
  if (distribution %in% c("normal", "uniform", "beta")) {
    result <- num.commas == 1
  }
  if (distribution %in% c("exponential", "poisson")) {
    result <- num.commas == 0
  }
  result
}


# Assign Stan distribution integer ----------------------------------------
#' Return the integer value for STAN that corresponds to each supported distribution
#'
#' Currently supported distributions: normal (0), uniform (1), exponential (2).
#'
#' Used internally, in \code{\link{prep_init_list}}.
#'
#' @keywords internal
#'
#' @param distribution distribution
#'
#' @return Integer value corresponding to chosen distribution
#'
assign_stan_dist_int <- function(distribution) {
  if (distribution == "normal") {
    result <- as.integer(0)
  }
  if (distribution == "uniform") {
    result <- as.integer(1)
  }
  if (distribution == "exponential") {
    result <- as.integer(2)
  }
  result
}

