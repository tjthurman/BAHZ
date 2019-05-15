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


plot_cline <- function(stanfit) {

  assertthat::assert_that(class(stanfit)[1] == "stanfit",
                          msg = "Model fit to plot from must be of class stanfit")

  # Get summary of the model
  summary <- bahz::cline_summary(stanfit, show.all = T)

  # Figure out if increasing or decreasing
  ps <- summary %>%
    filter(str_detect(param, "^p\\["))
  if (ps$mean[1] < ps$mean[dim(ps)[1]]) {
    decreasing <- F
  } else {
    decreasing = T
  }

  # Extract the cline parameters used
  cline_params <- summary %>%
    filter(param %in% c("center", "width",
                        "pmin", "pmax",
                        "deltaL", "deltaR", "deltaM",
                        "tauL", "tauR", "tauM")) %>%
    select(param, mean)

  # Get params to pass to general_cline_equation
  center <- ifelse("center" %in% cline_params$param,
                   cline_params$mean[which(cline_params$param == "center")],
                   NULL)
  width <- ifelse("width" %in% cline_params$param,
                   cline_params$mean[which(cline_params$param == "width")],
                   NULL)
  pmin <- ifelse("pmin" %in% cline_params$param,
                  cline_params$mean[which(cline_params$param == "pmin")],
                  NULL)
  pmax <- ifelse("pmax" %in% cline_params$param,
                  cline_params$mean[which(cline_params$param == "pmax")],
                  NULL)
  if ("deltaM" %in% cline_params$param) {
    deltaL <- cline_params$mean[which(cline_params$param == "deltaM")]
    deltaR <- cline_params$mean[which(cline_params$param == "deltaM")]
  } else if ("deltaL" %in% cline_params$param) {
    deltaL <- cline_params$mean[which(cline_params$param == "deltaL")]
    deltaR <- cline_params$mean[which(cline_params$param == "deltaR")]
  } else {
    deltaL <- NULL
    deltaR <- NULL
  }
  if ("tauM" %in% cline_params$param) {
    tauL <- cline_params$mean[which(cline_params$param == "tauM")]
    tauR <- cline_params$mean[which(cline_params$param == "tauM")]
  } else if ("tauL" %in% cline_params$param) {
    tauL <- cline_params$mean[which(cline_params$param == "tauL")]
    tauR <- cline_params$mean[which(cline_params$param == "tauR")]
  } else {
    tauL <- NULL
    tauR <- NULL
  }

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

  tibble(transectDist = xrange,
         p = y)

}
