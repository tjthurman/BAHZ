#' Calculate allele frequency at a position on a cline
#'
#' A flexible function to evaluate hybrid zone genetic cline equations, with or
#' without introgression tails.
#'
#' DETAILS TO BE ADDED. SPECIFICALLY, SOME NOTES ON EQUATIONS.
#'
#' @param transectDist The position at which to evaluate the cline equation.
#' @param decrease Is the cline decreasing in frequency? \code{TRUE} or
#'   \code{FALSE}.
#' @param center The location of the cline center, in the same distance units as
#'   \code{transectDist}. Numeric, must be greater than 0.
#' @param width The width of the cline, in the same distance units as
#'   \code{transectDist}. Numeric, must be greater than 0.
#' @param pmin,pmax Optional. The minimum and maximum allele frequency values
#'   in the tails of the cline. Default values are \code{0} and \code{1}, respectively.
#'   Must be between 0 and 1 (inclusive). Numeric.
#' @param deltaL,tauL Optional delta and tau parameters which describe the left
#'   exponential tail. Must supply both to generate a tail. Default is \code{NULL} (no
#'   tails). Numeric. tauL must be between 0 and 1 (inclusive).
#' @param deltaR,tauR Optional delta and tau parameters which describe the right
#'   exponential tail. Must supply both to generate a tail. Default is \code{NULL} (no
#'   tails). Numeric. tauR must be between 0 and 1 (inclusive).
#'
#' @return The allele frequency at that point of the cline.
#'
#' @export
#'
#' @examples
#' # Calculate the allele frequency at x = 100 for an increasing cline
#' with a center at x = 125, a width of 40, and no tails.
#'
#' general_cline_eqn(transectDist = 100, decrease = F,
#'                   center = 125, width = 40)
#'
#' # Calculate the allele frequency at x = 84 for a decreasing cline
#' with a center at 76, a width of 22, and a right tail
#' with deltaR = 7.5, tauR = 0.45.
#'
#' general_cline_eqn(transectDist = 84, decrease = F,
#'                  center = 76, width = 22,
#'                  deltaR = 7.5, tauR = .45)
#'

general_cline_eqn <- function(transectDist, decrease,
                              center, width,
                              pmin = 0, pmax = 1,
                              deltaL = NULL, tauL = NULL,
                              deltaR = NULL, tauR = NULL) {

  # Start with an ungodly amount of argument checking
  # check transect distance is a numeric vector
  assertthat::assert_that(is.vector(transectDist) == T, msg = "transect_distances must be a vector")
  assertthat::assert_that(is.numeric(transectDist) == T, msg = "transect_distances must be of type numeric")

  # Check decrease is T/F
  assertthat::assert_that(is.logical(decrease) == T, msg = "decrease must be either TRUE (T) or FALSE (F)")

  # Check the cline parameters are numeric vectors of length 1
  for (num.arg in alist(center, width, pmin, pmax, deltaL, deltaR, tauL, tauR)) {
    if (is.null(eval(num.arg)) == F) {
      assertthat::assert_that(is.vector(eval(num.arg)) == T, msg = paste(num.arg, "must be a vector", sep = " "))
      assertthat::assert_that(is.numeric(eval(num.arg)) == T, msg = paste(num.arg, "must be numeric", sep = " "))
      assertthat::assert_that(length(eval(num.arg)) == 1, msg = paste(num.arg, "must be of length 1", sep = " "))
    }
  }
  # Center and width must be greater than 0
  assertthat::assert_that(center >= 0, msg = "center must be greater than 0")
  assertthat::assert_that(width >= 0, msg = "width must be greater than 0")

  # Pmin, pmax, tauL, and tauR must be between 0 and 1
  for (num.arg in alist(pmin, pmax, tauL, tauR)) {
    if (is.null(eval(num.arg)) == F) {
    assertthat::assert_that(eval(num.arg) >= 0, msg = paste(num.arg, " must be between 0 and 1 (inclusive)", sep = ""))
    assertthat::assert_that(eval(num.arg) <= 1, msg = paste(num.arg, " must be between 0 and 1 (inclusive)", sep = ""))
    }
  }

  # Check to make sure both delta and tau are supplied for a given set of tails
  assertthat::assert_that(sum(is.null(deltaL), is.null(tauL)) %in% c(0,2),
              msg = "If using deltaL or tauL, must supply both of them")
  assertthat::assert_that(sum(is.null(deltaR), is.null(tauR)) %in% c(0,2),
              msg = "If using deltaR or tauR, must supply both of them")

  # The cline equations are written for increasing clines
  # With two alleles, a decreasing cline is simply
  # 1- increasing
  # Will do that at the end, after computng p.

  # Calculate p, using the proper equation based on the parameters prodived
  if (sum(is.null(deltaL), is.null(deltaR)) == 2) { # if no tails
    prop <- (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
  }
  else if (sum(is.null(deltaL), is.null(deltaR)) == 0) {# If both tails
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
  else if (is.null(deltaL) == F) { # If left tail only
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


  if (decrease == T) {
    p <- pmin + (pmax - pmin)*(1 - prop)
  }
  if (decrease == F) {
    p <- pmin + (pmax - pmin)*prop
  }
  p
}
