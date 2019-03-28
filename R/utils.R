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
#' A function that corrects Fis estimates which are negative or NaN, turning them
#' to 0. Reads in a data frame with an Fis column. Meant to be part of a dplyr
#' pipeline.
#'
#' Used internally, in \code{\link{sim_geno_cline}}.
#'
#' @keywords internal
#'
#' @param .df A data frame, containing a column named Fis.
#'
#' @return The supplied data frame, with a corrected Fis column.
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
#' @param prior_file Filepath to the prior file.
#'
#' @return A named list containing the specifed priors.
#'
#' @examples
#' \dontrun{
#' parse_prior_file(path/to/priors.yaml)
#' }
#'

parse_prior_file <- function(prior_file) {
  path_to_prior <-
    file.path(normalizePath(prior_file), fsep = .Platform$file.sep)
  priors <- yaml::yaml.load_file(path_to_prior, as.named.list = T)

  # Move this check up to the prior parsing? Will be more useful there, can say which
  # one is causing the issue
  assertthat::assert_that(is.character(unlist(priors)),
                          msg = "Problem parsing priors. Check your prior file!")

  # Check for proper length and that order is correct
  if (length(names(priors)) != 11) {
    stop("Incorrect number of priors config file! Double check your file against the template!")
  }
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
#' @description Internal functions used in \code{\link{prep_init_list}}
#' and \code{\link{prep_prior_list}} to extract
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
#' @param string The string to extract a value from.
#'
#' @return A numeric value.
#'
#' @examples
#' \dontrun{
#' extract_first("uniform(32,45)") # returns 32
#' extract_last("uniform(32,45)") # returns 45
#' extract_only("exponential(32)") # returns 32
#' }

NULL


#'
#' Extract values for prior distributions
#'
#' @rdname extractValue
#'
#'

extract_first <- function(string) {

  chr.res <- stringr::str_extract(string, "\\(.*,") %>%
    stringr::str_remove_all("[(,)]") %>%
    stringr::str_replace_all(pattern = " ", "")

  assertthat::assert_that(stringr::str_count(chr.res, "\\.") <= 1,
                          msg = "Could not parse prior, too many decimal places. Check your prior file!")
  if (suppressWarnings(is.na(as.numeric(chr.res)))) {
    stop("Could not coerce prior to a numeric value. Check your prior file!")
  }
  res <- as.numeric(chr.res)
  assertthat::assert_that(length(res) == 1, msg = "Could not parse prior. Check your prior file!")
  res
}


#' @rdname extractValue
#'
#'

extract_last <- function(string) {
  chr.res <- stringr::str_extract(string, ",.*\\)") %>%
    stringr::str_remove_all("[(,)]") %>%
    stringr::str_replace_all(pattern = " ", "")
  assertthat::assert_that(stringr::str_count(chr.res, "\\.") <= 1,
                          msg = "Could not parse prior, too many decimal places. Check your prior file!")
  if (suppressWarnings(is.na(as.numeric(chr.res)))) {
    stop("Could not coerce prior to a numeric value. Check your prior file!")
  }
  res <- as.numeric(chr.res)
  assertthat::assert_that(length(res) == 1, msg = "Could not parse prior. Check your prior file!")
  res
}


#'
#' Extract values for prior distributions
#'
#' @rdname extractValue
#'
#'
extract_only <- function(string) {
  chr.res <- stringr::str_extract(string, "\\(.*\\)") %>%
    stringr::str_remove_all("[()]") %>%
    stringr::str_replace_all(pattern = " ", "")
  assertthat::assert_that(stringr::str_count(chr.res, ",") == 0,
                          msg = "Could not parse prior, comma present in a distribution that requires a single value. Check your prior file!")
  assertthat::assert_that(stringr::str_count(chr.res, "\\.") <= 1,
                          msg = "Could not parse prior, too many decimal places. Check your prior file!")
  if (suppressWarnings(is.na(as.numeric(chr.res)))) {
    stop("Could not coerce prior to a numeric value. Check your prior file!")
  }
  res <- as.numeric(chr.res)
  assertthat::assert_that(length(res) == 1, msg = "Could not parse prior. Check your prior file!")
  res
}


# Check that priors are supported ---------------------------------------

#' Check that the distribution specified for a given parameter is supported
#'
#' Used internally, in \code{\link{prep_prior_list}}. Checks that the prior
#' distribution specified for a particular cline parameter is currently
#' supported.
#'
#' @keywords internal
#'
#' @param parameter The parameter for which a prior distribution is specified,
#'   as a character string.
#'
#' @param distribution The prior distribution for that parameter, as a character
#'   string.
#'
#' @return TRUE/FALSE: is the distribution supported for that parameter?
#'
#' @examples
#' \dontrun{
#' check_prior_supported("center", "normal") # returns T
#' check_prior_supported("center", "beta") # returns F
#' check_prior_supported("deltaL", "normal") # returns F
#' }


check_prior_supported <- function(parameter, distribution) {
  result <- FALSE # "Failsafe" to false
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
#' @param distribution the prior distribution, as a character string.
#'
#' @param string the string specifiying the parameter values, parsed from the prior_config file.
#'
#' @return TRUE/FALSE: is the prior specified correctly?
#'
#' @examples
#' \dontrun{
#' check_prior_specification("normal", "normal(32,45)") # returns T
#' check_prior_specification("uniform", "uniform(32)") # returns F
#' check_prior_specification("exponential", "exponential(0.5)") # returns T
#' }
#'

check_prior_specification <- function(distribution, string) {
  result <- F # "Failsafe" to false
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
#' @param distribution The specified distribution, as a character string.
#'
#' @return Integer value corresponding to chosen distribution.
#'
#' @examples
#' \dontrun{
#' assign_stan_dist_int("normal") # returns 0
#' assign_stan_dist_int("uniform") # returns 1
#' assign_stan_dist_int("exponential") # returns 2
#' }
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


# Check init chain --------------------------------------------------------
#' Check whether the random initial values for a chain are appropriate.
#'
#' Checks the initial values before they are passed to Stan to make sure they are appropriate.
#'
#' Used internally, in \code{\link{prep_init_list}}.
#'
#' Checks that: width is positive, delta parameters are positive, tau
#' parameters, pmin/pmax, and f are between 0 and 1, and that pmin is less than
#' pmax.
#'
#' @keywords internal
#'
#' @param single.init.list The list of initial values to be checked.
#'
#' @return NULL, if no problems, otherwise a vector of parameters which have
#'   inappropriate initial values.
#'
#' @examples
#' \dontrun{
#' assign_stan_dist_int("normal") # returns 0
#' assign_stan_dist_int("uniform") # returns 1
#' assign_stan_dist_int("exponential") # returns 2
#' }
#'
check_init_chain <- function(single.init.list) {
  problems <- NULL
  for (j in 1:length(single.init.list)) {
    if (names(single.init.list)[j] == "width") {
      if (single.init.list[[j]] < 0) {
        problems <- c(problems, "width")
      }
    }
    if (names(single.init.list)[j] == "deltaL") {
      if (single.init.list[[j]] < 0) {
        problems <- c(problems, "deltaL")
      }
    }
    if (names(single.init.list)[j] == "deltaR") {
      if (single.init.list[[j]] < 0) {
        problems <- c(problems, "deltaR")
      }
    }
    if (names(single.init.list)[j] == "deltaM") {
      if (single.init.list[[j]] < 0) {
        problems <- c(problems, "deltaM")
      }
    }
    if (names(single.init.list)[j] == "tauL") {
      if (dplyr::between(single.init.list[[j]], 0, 1) == F) {
        problems <- c(problems, "tauL")
      }
    }
    if (names(single.init.list)[j] == "tauR") {
      if (dplyr::between(single.init.list[[j]], 0, 1) == F) {
        problems <- c(problems, "tauR")
      }
    }
    if (names(single.init.list)[j] == "tauM") {
      if (dplyr::between(single.init.list[[j]], 0, 1) == F) {
        problems <- c(problems, "tauM")
      }
    }
    if (names(single.init.list)[j] == "pmin") {
      if (dplyr::between(single.init.list[[j]], 0, 1) == F) {
        problems <- c(problems, "pmin")
      }
    }
    if (names(single.init.list)[j] == "pmax") {
      if (dplyr::between(single.init.list[[j]], 0, 1) == F) {
        problems <- c(problems, "pmax")
      }
    }
    if (names(single.init.list)[j] == "f") {
      if (dplyr::between(single.init.list[[j]], 0, 1) == F) {
        problems <- c(problems, "f")
      }
    }
    if (names(single.init.list)[j] == "pmax") {
      if (single.init.list$pmin > single.init.list[[j]]) {
        problems <- c(problems, "pmin", "pmax")
      }
    }
  }
  problems
}



