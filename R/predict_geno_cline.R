#' Pot thing
#'
#'
#' @importClassesFrom rstan stanfit
#'
#' @param stanfit A \code{\linkS4class{stanfit}} object holding your model results.
#'
#' @param data The dataframe with your cline data.
#'
#' @param num.out Optional, the number of x values (distances) at which to
#'   evaluate the cline. By default, does twice the length of the cline, but you
#'   can specify more or fewer. Too few may leade to a jagged-looking cline.
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' predict_genocline(yourStanfit, data)
#'
#' # If you want to specify only 100 data points for plotting
#' predict_genocline(yourStanfit, data, num.out = 100)
#' }


predict_geno_cline <- function(stanfit, data, num.out = NULL) {

  # Check arguments
  # add checks for data
  # add checks for num.out
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


  # Figure out starting and ending x values
  # Give 2% visual padding

  start.x <- as.integer(min(data$transectDist))
  end.x <- as.integer(max(data$transectDist))
  range <- end.x - start.x
  start.pad <- as.integer(start.x - (0.02*(range)))
  end.pad <- as.integer(end.x + (0.02*(range)))

  xrange <- seq(from = start.pad, to = end.pad,
                length.out = ifelse(is.null(num.out) == T,
                                    range*2,
                                    num.out))

  #Pass those in to a loop with the general_cline_equation
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


  dplyr::tibble(transectDist = xrange,
         p = y)
}
