#' Predict allele frequencies for a genetic cline
#'
#' Uses the cline parameters from your stanfit object (the mean of the posterior
#' distribtution of each parameter) to predict the allele frequency at the
#' distance(s) along the transect specified by the user.
#'
#' @importClassesFrom rstan stanfit
#'
#' @param stanfit A \code{\linkS4class{stanfit}} object holding your model results.
#'
#' @param distance The x value(s) (distance along the transect) at which to
#'   predict allele frequencies for the fitted cline. Must be a numeric vector.
#'
#' @return A data frame with two columns:
#'  \itemize{
#'  \item transectDist: the distances along the transect at which the cline
#'  equation was evaluated. Will be the same values as were supplied in the
#'  distance argument. Numeric.
#'  \item p: the predict allele frequency at each distance, according to the
#'  cline equation. Numeric.
#'  }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Predict at one distance
#' predict_geno_cline(yourStanfit, distance = 50)
#'
#' # Predict at a range of distances
#' predict_geno_cline(yourStanfit, distance = 0:100)
#'
#' }


predict_geno_cline <- function(stanfit, distance) {

  # Check arguments
  assertthat::assert_that(class(stanfit)[1] == "stanfit",
                          msg = "Model object from which to plot must be of class stanfit")
  assertthat::assert_that(is.vector(distance),
                          msg = paste("distance must be a vector"))
  assertthat::assert_that(is.numeric(distance),
                          msg = paste("distance must be of type numeric, not ",
                                      typeof(distance), sep = ""))

  # Get summary of the model
  summ <- bahz::cline_summary(stanfit, show.all = T)

  # Figure out if increasing or decreasing
  ps <- summ %>%
    dplyr::filter(stringr::str_detect(.$param, "^p\\["))
  if (ps$mean[1] < ps$mean[dim(ps)[1]]) {
    decreasing <- F
  } else {
    decreasing = T
  }

  # Get params to pass to general_cline_equation
  # First 4 are always there
  center <- ifelse("center" %in% summ$param,
                   summ$mean[which(summ$param == "center")],
                   NULL)
  width <- ifelse("width" %in% summ$param,
                   summ$mean[which(summ$param == "width")],
                   NULL)
  pmin <- ifelse("pmin" %in% summ$param,
                  summ$mean[which(summ$param == "pmin")],
                  NULL)
  pmax <- ifelse("pmax" %in% summ$param,
                  summ$mean[which(summ$param == "pmax")],
                  NULL)
  # Next few are optional based on tails
  if ("deltaM" %in% summ$param) { # If mirror
    deltaL <- summ$mean[which(summ$param == "deltaM")]
    deltaR <- summ$mean[which(summ$param == "deltaM")]
    tauL <- summ$mean[which(summ$param == "tauM")]
    tauR <- summ$mean[which(summ$param == "tauM")]
  } else if ("deltaL" %in% summ$param) { # if there's a left
    deltaL <- summ$mean[which(summ$param == "deltaL")]
    tauL <- summ$mean[which(summ$param == "tauL")]
    if ("deltaR" %in% summ$param) { # if there's also a right, independent
      deltaR <- summ$mean[which(summ$param == "deltaR")]
      tauR <- summ$mean[which(summ$param == "tauR")]
    } else { # else left only
      deltaR <- NULL
      tauR <- NULL
    }
  } else if ("deltaR" %in% summ$param) { # if right only
      deltaR <- summ$mean[which(summ$param == "deltaR")]
      tauR <- summ$mean[which(summ$param == "tauR")]
      deltaL <- NULL
      tauL <- NULL
  } else { # no tails
    deltaL <- NULL
    deltaR <- NULL
    tauL <- NULL
    tauR <- NULL
  }


  #Pass those to general cline equation.
  #in sapply for now, will vectorize that soon.
  y <- sapply(distance, FUN = bahz::general_cline_eqn,
              decrease = decreasing,
              center = center,
              width = width,
              pmin = pmin,
              pmax = pmax,
              deltaL = deltaL,
              tauL = tauL,
              deltaR = deltaR,
              tauR = tauR,
              USE.NAMES = F)


  # Some ideas on how to add lines for every draw from the posterior.
  # But, pretty slow. Would want to find a way to vectorize/apply over
  # the posterior matrix.

  # posterior <- as.matrix(stanfit)
  # y_post <- matrix(NA, nrow = length(xrange), ncol = dim(posterior)[1])
  # for (sample in 1:dim(posterior)[1]) {
  #   i <- 1
  #   for (x in xrange) {
  #     y_post[i, sample] <- bahz::general_cline_eqn(transectDist = x, decrease = decreasing,
  #                                     center = posterior[sample, which(colnames(posterior) == "center")],
  #                                     width = posterior[sample, which(colnames(posterior) == "width")],
  #                                     pmin = posterior[sample, which(colnames(posterior) == "pmin")],
  #                                     pmax = posterior[sample, which(colnames(posterior) == "pmax")],
  #                                     deltaL = NULL,
  #                                     tauL = NULL,
  #                                     deltaR = NULL,
  #                                     tauR = NULL)
  #     i <- i + 1
  #   }
  # }


  data.frame(transectDist = distance,
         p = y, stringsAsFactors = F)
}
