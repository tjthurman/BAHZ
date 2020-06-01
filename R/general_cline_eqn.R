#' Calculate allele frequencies for a cline
#'
#' A flexible function to evaluate hybrid zone genetic cline equations, with or
#' without introgression tails.
#'
#' This function evaluates hybrid zone cline equations of the general form:
#'
#'  \deqn{p = (e^(4*((x-center))/(width))/((1+ e^(4*((x-center))/(wwidth)))}
#'
#'  Where \eqn{center} is the center of the cline, \eqn{width} is the width of
#'  the cline, \eqn{x} is the distance along the transect, and \eqn{p} is the
#'  expected frequency of the allele.
#'
#'  This function can also include introgression tails by specifying delta and
#'  tau paramters. For the left introgression tail, when \eqn{x \le (center -
#'  \delta)}:
#'
#'  \deqn{p = 1/(1 + e^(4(\delta)/(w)))*\exp((4\tau(x-c+\delta)/width)/(1+e^(-4(\delta)/(w))))}
#'
#'  For the right introgression tail, when \eqn{x \ge (center + \delta)}:
#'
#'  \deqn{p = 1- 1/(1 + e^(4(\delta)/(w)))*\exp((-4\tau(x-c-\delta)/width)/(1+e^(-4(\delta)/(w))))}
#'
#'  The above equations all describe an increasing cline with minimum and
#'  maximum values of 0 and 1. This equation can be rescaled for other
#'  minima and maxima:
#'
#'  \deqn{p' = p_{min} + (p_{max} - p_{min})*p}
#'
#'  And for clines which decrease in frequency:
#'
#'  \deqn{p' = p_{min} + (p_{max} - p_{min})*(1-p)}
#'
#'  The function will automatically use the appropriate equation based on
#'  arguments supplied by the user.
#'
#'
#' @param transectDist The position(s) at which to evaluate the cline equation.
#'   Must be a numeric vector.
#' @param decrease Is the cline decreasing in frequency? \code{TRUE} or
#'   \code{FALSE}.
#' @param center The location of the cline center, in the same distance units as
#'   \code{transectDist}. Numeric.
#' @param width The width of the cline, in the same distance units as
#'   \code{transectDist}. Numeric, must be greater than 0.
#' @param pmin,pmax Optional. The minimum and maximum allele frequency values in
#'   the tails of the cline. Default values are \code{0} and \code{1},
#'   respectively. Must be between 0 and 1 (inclusive). Numeric.
#' @param deltaL,tauL Optional delta and tau parameters which describe the left
#'   exponential introgression tail. Must supply both to generate a tail. Default is
#'   \code{NA} (no tails). Numeric. tauL must be between 0 and 1 (inclusive).
#' @param deltaR,tauR Optional delta and tau parameters which describe the right
#'   exponential introgression tail. Must supply both to generate a tail. Default is
#'   \code{NA} (no tails). Numeric. tauR must be between 0 and 1 (inclusive).
#'
#' @return The result of evaluating a cline equation with the specified
#'   parameters at the specified distance(s) (that is, \eqn{p'} in the notation of the
#'   equations above). A numeric vector equal in length to the transectDist argument.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate the allele frequency at x = 100 for an increasing cline
#' # with a center at x = 125, a width of 40, and no tails.
#'
#' general_cline_eqn(transectDist = 100, decrease = F,
#'                   center = 125, width = 40)
#'
#' # Calculate the allele frequency at x = 84 and x = 100 for a decreasing cline
#' # with a center at 76, a width of 22, and a right tail
#' # with deltaR = 7.5, tauR = 0.45.
#'
#' general_cline_eqn(transectDist = c(84,100), decrease = TRUE,
#'                  center = 76, width = 22,
#'                  deltaR = 7.5, tauR = .45)
#' }
#'
#'


general_cline_eqn <- function(transectDist, decrease = c(TRUE, FALSE),
                                      center, width,
                                      pmin = 0, pmax = 1,
                                      deltaL = NA, tauL = NA,
                                      deltaR = NA, tauR = NA) {
  # Start with an ungodly amount of argument checking
  # check transect distance is a numeric vector
  assertthat::assert_that(is.vector(transectDist) == T, msg = "transect_distances must be a vector")
  assertthat::assert_that(is.numeric(transectDist) == T, msg = "transect_distances must be of type numeric")

  # Check decrease is T/F
  assertthat::assert_that(is.logical(decrease) == T, msg = "decrease must be either TRUE (T) or FALSE (F)")

  # Check the cline parameters are of length 1
  for (num.arg in alist(center, width, pmin, pmax, deltaL, deltaR, tauL, tauR)) {
    if (is.null(eval(num.arg)) == F) {
      assertthat::assert_that(length(eval(num.arg)) == 1, msg = paste(num.arg, "must be of length 1", sep = " "))
    }
  }
  # Check that the cline parameters are numeric and vectors
  for (num.arg in alist(center, width, pmin, pmax, deltaL, deltaR, tauL, tauR)) {
    if (is.na(eval(num.arg)) == F) {
      assertthat::assert_that(is.vector(eval(num.arg)) == T, msg = paste(num.arg, "must be a vector", sep = " "))
      assertthat::assert_that(is.numeric(eval(num.arg)) == T, msg = paste(num.arg, "must be numeric", sep = " "))
    }
  }
  # Width must be greater than 0
  assertthat::assert_that(width >= 0, msg = "width must be greater than 0")

  # Pmin, pmax, tauL, and tauR must be between 0 and 1
  for (num.arg in alist(tauL, tauR)) {
    if (is.na(eval(num.arg)) == F) {
      assertthat::assert_that(eval(num.arg) >= 0, msg = paste(num.arg, " must be between 0 and 1 (inclusive)", sep = ""))
      assertthat::assert_that(eval(num.arg) <= 1, msg = paste(num.arg, " must be between 0 and 1 (inclusive)", sep = ""))
    }
  }

  # Check to make sure both delta and tau are supplied for a given set of tails
  assertthat::assert_that(sum(is.na(deltaL), is.na(tauL)) %in% c(0,2),
                          msg = "If using deltaL or tauL, must supply both of them")
  assertthat::assert_that(sum(is.na(deltaR), is.na(tauR)) %in% c(0,2),
                          msg = "If using deltaR or tauR, must supply both of them")

  # Then call internal_cline_eqn within sapply

  sapply(transectDist, FUN = internal_cline_eqn,
         decrease = decrease,
         center = center,
         width = width,
         pmin = pmin,
         pmax = pmax,
         deltaL = deltaL,
         tauL = tauL,
         deltaR = deltaR,
         tauR = tauR,
         USE.NAMES = F)
}


#' Calculate allele frequencies for a single point along a cline
#'
#' A flexible function to evaluate hybrid zone genetic cline equations at a
#' single point. Used internally by \code{\link{general_cline_eqn}}, which
#' calculates results for multiple points. See that help page for all the math
#' details.
#'
#' @param transectDist The position at which to evaluate the cline equation.
#'   Must be a numeric vector.
#' @param decrease Is the cline decreasing in frequency? \code{TRUE} or
#'   \code{FALSE}.
#' @param center The location of the cline center, in the same distance units as
#'   \code{transectDist}. Numeric.
#' @param width The width of the cline, in the same distance units as
#'   \code{transectDist}. Numeric, must be greater than 0.
#' @param pmin,pmax Optional. The minimum and maximum allele frequency values in
#'   the tails of the cline. Default values are \code{0} and \code{1},
#'   respectively. Must be between 0 and 1 (inclusive). Numeric.
#' @param deltaL,tauL Optional delta and tau parameters which describe the left
#'   exponential introgression tail. Must supply both to generate a tail. Default is
#'   \code{NA} (no tails). Numeric. tauL must be between 0 and 1 (inclusive).
#' @param deltaR,tauR Optional delta and tau parameters which describe the right
#'   exponential introgression tail. Must supply both to generate a tail. Default is
#'   \code{NA} (no tails). Numeric. tauR must be between 0 and 1 (inclusive).
#'
#' @keywords internal
#'
#' @return The result of evaluating a cline equation with the specified
#'   parameters at the specified distance (that is, \eqn{p'} in the notation of the
#'   equations above). A numeric vector of length 1.


internal_cline_eqn <- function(transectDist, decrease = c(TRUE, FALSE),
                              center, width,
                              pmin = 0, pmax = 1,
                              deltaL = NA, tauL = NA,
                              deltaR = NA, tauR = NA) {



  # The cline equations are written for increasing clines
  # With two alleles, a decreasing cline is simply
  # 1- increasing
  # Will do that at the end, after computng p.

  # Calculate p, using the proper equation based on the parameters prodived
  if (sum(is.na(deltaL), is.na(deltaR)) == 2) { # if no tails
    prop <- (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
  }
  else if (sum(is.na(deltaL), is.na(deltaR)) == 0) {# If both tails
    if (transectDist <= center - deltaL) { # and we're in the left tail
      # use the left tail equation
      prop <- (1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist - center + deltaL)/width)/(1 + exp(-4*deltaL/width)))
    }
    else if (transectDist >= center + deltaR) { # and we're in the right tail
      # use the right tail equation
      prop <- (1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist - center - deltaR)/width)/(1 + exp(-4*deltaR/width))))
    }
    else {# Otherwise, we're in the sigmoid center
      prop <- (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
    }
  }
  else if (is.na(deltaL) == F) { # If left tail only
    if (transectDist <= center - deltaL) { # and we're in the left tail
      # use the left tail equation
      prop <- (1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist - center + deltaL)/width)/(1 + exp(-4*deltaL/width)))
    }
    else {# Otherwise, we're in the sigmoid center
      prop <- (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
    }
  }
  else { # must be right tail cline
    if(transectDist >= center + deltaR) { # and we're in the right tail
      # use the right tail equation
      prop <- (1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist - center - deltaR)/width)/(1 + exp(-4*deltaR/width))))
    }
    else {# Otherwise, we're in the sigmoid center or there are no tails
      prop <- (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
    }
  }

  # The exponential function in R will got to Inf once it hits 710.
  # So, if (distance-center)/width gets to be 177.5, the exp() function goes to Inf
  # This turns prop into NaN, when it should be 1.
  # Will just add a check for this.
  if (is.nan(prop) == T) {
    prop <- 1
  }

  # Rescale to pmin/pmax and do increase/decrease
  if (decrease == T) {
    p <- pmin + (pmax - pmin)*(1 - prop)
  }
  if (decrease == F) {
    p <- pmin + (pmax - pmin)*prop
  }
  p
}
