#' Pot thing
#'
#'
#' @importClassesFrom rstan stanfit
#'
#' @param stanfit A \code{\linkS4class{stanfit}} object holding your model results.
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_cline(yourStanfit)
#' }


plot_cline <- function(stanfit) {

  assertthat::assert_that(class(stanfit)[1] == "stanfit",
                          msg = "Model object from which to plot must be of class stanfit")

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

  #Pass those in to a loop with the general_cline_equation
  xrange <- -300:300
  y <- NULL
  i <- 1
  for (x in xrange) {
    y[i] <- bahz::general_cline_eqn(transectDist = x, decrease = decreasing,
                              center = center,
                              width = width,
                              pmin = pmin,
                              pmax = pmax,
                              deltaL = deltaL,
                              tauL = tauL,
                              deltaR = deltaR,
                              tauR = tauR)
    i <- i + 1
  }

  dplyr::tibble(transectDist = xrange,
         p = y)

}
