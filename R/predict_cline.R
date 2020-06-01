#' Predict allele frequencies or mean phenotypes for a cline
#'
#' Uses the cline parameters from your stanfit object to predict the
#' allele frequency (for genetic clines) or the mean phenotype value
#' (for phenotypic clines) at the distance(s) along the transect specified by the user.
#' Returns best-fit cline values estimated using both the mean of the
#' posterior distribution of each parameter and the median of the posterior
#' distribution of each parametet. Optionally, can also calculate credible
#' intervals around these best-fit estimates. To save computation time, this function use the
#' \code{\link[memoise]{memoise}} package to save past results and avoid recalculations. To clear the cache
#' of saved results, use clear.cache = F or \code{\link{clear_bahz_cache}}.
#'
#' @importClassesFrom rstan stanfit
#'
#' @import progress
#'
#' @import memoise
#'
#' @param stanfit A \code{\linkS4class{stanfit}} object holding your model results.
#'
#' @param distance The x value(s) (distance along the transect) at which to
#'   predict allele frequencies or mean phenotypes for the fitted cline. Must be
#'   a numeric vector.
#'
#' @param confidence Calculate credible intervals around the cline? TRUE or FALSE, default FALSE.
#'
#' @param prob The probability interval to calculate for the cline. Default is .95. Numeric,
#'   between 0 and 1.
#'
#' @param method The method for calculating credible intervals. Either "ET" for
#'   equal-tail probability intervals, or "HPDI" for highest posterior density
#'   intervals. Default is "HPDI".
#'
#' @param progress Show progress bar when calculating credible intervals? TRUE or FALSE, default TRUE.
#'
#' @param clear.cache Clear the cache of saved results to ensure recalculation of predicted cline?
#' TRUE or FALSE, default FALSE.
#'
#' @return A data frame with either 2 (confidence = F) or 4 (confidence = T) columns:
#'  \itemize{
#'  \item transectDist: the distances along the transect at which the cline
#'  equation was evaluated. Will be the same values as were supplied in the
#'  distance argument. Numeric.
#'  \item p_mean: the predicted allele frequency or mean phenotype at each distance, using the
#'  posterior mean for each parameter of the cline equation. Numeric.
#'  \item p_median: the predicted allele frequency or mean phenotype at each distance, using the
#'  posterior median for each parameter of the cline equation. Numeric.
#'  \item lower: Only if confidence = T. The lower limit of the credible interval. The column name
#'  will include the probability value used for calculating the CI.
#'  \item upper: Only if confidence = T. The upper limit of the credible interval. The column name
#'  will include the probability value used for calculating the CI.
#'  }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Predict at one distance
#' predict_cline(yourStanfit, distance = 50)
#'
#' # Predict at a range of distances
#' predict_cline(yourStanfit, distance = 0:100)
#'
#' }

predict_cline <- memoise::memoise(function(stanfit,
                                           distance,
                                           confidence = F,
                                           prob = 0.95,
                                           method = "HPDI",
                                           progress = T,
                                           clear.cache = F) {

  # Check arguments
  assertthat::assert_that(is.vector(distance),
                          msg = paste("distance must be a vector"))
  assertthat::assert_that(is.numeric(distance),
                          msg = paste("distance must be of type numeric, not ",
                                      typeof(distance), sep = ""))
  assertthat::assert_that(is.numeric(prob) == T, msg = "prob must be numeric")
  assertthat::assert_that(length(prob) == 1, msg = "prob must be of length 1")
  assertthat::assert_that(prob <= 1, msg = "prob must be between 0 and 1")
  assertthat::assert_that(prob > 0, msg = "prob must be between 0 and 1")
  assertthat::assert_that((method %in% c("HPDI", "ET")) == T,
                          msg = "method must be either 'HPDI' or 'ET")
  assertthat::assert_that(is.logical(confidence) == T, msg = "confidence must be either TRUE or FALSE")
  assertthat::assert_that(is.logical(progress) == T, msg = "progress must be either TRUE or FALSE")
  assertthat::assert_that(is.logical(clear.cache) == T, msg = "clear.cache must be either TRUE or FALSE")

  if (clear.cache) {
    memoise::forget(bahz::predict_cline)
  }

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
  c.mean <- ifelse("center" %in% summ$param,
                   summ$mean[which(summ$param == "center")],
                   NULL)
  c.median <- ifelse("center" %in% summ$param,
                   summ$median[which(summ$param == "center")],
                   NULL)
  w.mean <- ifelse("width" %in% summ$param,
                   summ$mean[which(summ$param == "width")],
                   NULL)
  w.median <- ifelse("width" %in% summ$param,
                   summ$median[which(summ$param == "width")],
                   NULL)
  pmin.mean <- ifelse("pmin" %in% summ$param,
                  summ$mean[which(summ$param == "pmin")],
                  NULL)
  pmin.median <- ifelse("pmin" %in% summ$param,
                      summ$median[which(summ$param == "pmin")],
                      NULL)
  pmax.mean <- ifelse("pmax" %in% summ$param,
                  summ$mean[which(summ$param == "pmax")],
                  NULL)
  pmax.median <- ifelse("pmax" %in% summ$param,
                      summ$median[which(summ$param == "pmax")],
                      NULL)
  # Next few are optional based on tails
  if ("deltaM" %in% summ$param) { # If mirror
    dL.mean <- summ$mean[which(summ$param == "deltaM")]
    dL.median <- summ$median[which(summ$param == "deltaM")]

    dR.mean <- summ$mean[which(summ$param == "deltaM")]
    dR.median <- summ$median[which(summ$param == "deltaM")]

    tL.mean <- summ$mean[which(summ$param == "tauM")]
    tL.median <- summ$median[which(summ$param == "tauM")]

    tR.mean <- summ$mean[which(summ$param == "tauM")]
    tR.median <- summ$median[which(summ$param == "tauM")]

  } else if ("deltaL" %in% summ$param) { # if there's a left
    dL.mean <- summ$mean[which(summ$param == "deltaL")]
    dL.median <- summ$median[which(summ$param == "deltaL")]

    tL.mean <- summ$mean[which(summ$param == "tauL")]
    tL.median <- summ$median[which(summ$param == "tauL")]

    if ("deltaR" %in% summ$param) { # if there's also a right, independent
      dR.mean <- summ$mean[which(summ$param == "deltaR")]
      dR.median <- summ$median[which(summ$param == "deltaR")]

      tR.mean <- summ$mean[which(summ$param == "tauR")]
      tR.median <- summ$median[which(summ$param == "tauR")]

    } else { # else left only
      dR.mean <- NA
      dR.median <- NA
      tR.mean <- NA
      tR.median <- NA
    }
  } else if ("deltaR" %in% summ$param) { # if right only
    dR.mean <- summ$mean[which(summ$param == "deltaR")]
    dR.median <- summ$median[which(summ$param == "deltaR")]

    tR.mean <- summ$mean[which(summ$param == "tauR")]
    tR.median <- summ$median[which(summ$param == "tauR")]

    dL.mean <- NA
    dL.median <- NA
    tL.mean <- NA
    tL.median <- NA
  } else { # no tails
    dL.mean <- NA
    dL.median <- NA
    tL.mean <- NA
    tL.median <- NA
    dR.mean <- NA
    dR.median <- NA
    tR.mean <- NA
    tR.median <- NA
  }

  #Pass those to general cline equation to get the
  # best fit cline
  y_mean <- general_cline_eqn(transectDist = distance,
                         decrease = decreasing,
                         center = c.mean,
                         width = w.mean,
                         pmin = pmin.mean,
                         pmax = pmax.mean,
                         deltaL = dL.mean,
                         tauL = tL.mean,
                         deltaR = dR.mean,
                         tauR = tR.mean)
  y_median <- general_cline_eqn(transectDist = distance,
                              decrease = decreasing,
                              center = c.median,
                              width = w.median,
                              pmin = pmin.median,
                              pmax = pmax.median,
                              deltaL = dL.median,
                              tauL = tL.median,
                              deltaR = dR.median,
                              tauR = tR.median)



  # If we want to add confidence intervals.
  # Could possibly speed this up a bit by writing around the loop somehow.
  # Not sure.
  if (confidence == T) {
    # make col names for CI columns
    # tail <- (1 - prob) / 2
    low.name <- paste("low", prob, method, sep = "_")
    up.name <- paste("up", prob, method, sep = "_")

    posterior <- as.matrix(stanfit)
    y_post <- matrix(NA, nrow = length(distance), ncol = dim(posterior)[1])

    # Set up vectors of parameters for general cline_eqn.
    post_center <- posterior[ , which(colnames(posterior) == "center")]
    post_width <- posterior[ , which(colnames(posterior) == "width")]
    post_pmin <- posterior[ , which(colnames(posterior) == "pmin")]
    post_pmax <- posterior[ , which(colnames(posterior) == "pmax")]
    if ("deltaM" %in% colnames(posterior)) { # If mirror
      post_deltaL <- posterior[ , which(colnames(posterior) == "deltaM")]
      post_deltaR <- posterior[ , which(colnames(posterior) == "deltaM")]
      post_tauL <- posterior[ , which(colnames(posterior) == "tauM")]
      post_tauR <- posterior[ , which(colnames(posterior) == "tauM")]
    } else if ("deltaL" %in% colnames(posterior)) { # if there's a left
      post_deltaL <- posterior[ , which(colnames(posterior) == "deltaL")]
      post_tauL <- posterior[ , which(colnames(posterior) == "tauL")]
      if ("deltaR" %in% colnames(posterior)) { # if there's also a right, independent
        post_deltaR <-  posterior[ , which(colnames(posterior) == "deltaR")]
        post_tauR <- posterior[ , which(colnames(posterior) == "tauR")]
      } else { # else left only
        post_deltaR <- NA
        post_tauR <- NA
      }
    } else if ("deltaR" %in% colnames(posterior)) { # if right only
      post_deltaR <-  posterior[ , which(colnames(posterior) == "deltaR")]
      post_tauR <- posterior[ , which(colnames(posterior) == "tauR")]
      post_deltaL <- NA
      post_tauL <- NA
    } else { # no tails
      post_deltaL <- NA
      post_deltaR <- NA
      post_tauL <- NA
      post_tauR <- NA
    }
    if (progress) {
      pb <- progress::progress_bar$new(total = dim(posterior)[1],
                             format = " Posterior sample [:bar] :current/:total",
                             clear = F, width= 60)
      pb$tick(0)
    }


    for (sample in 1:dim(posterior)[1]) {
      y_post[ , sample] <- bahz::general_cline_eqn(transectDist = distance, decrease = decreasing,
                                                   center = post_center[sample],
                                                   width = post_width[sample],
                                                   pmin = post_pmin[sample],
                                                   pmax = post_pmax[sample],
                                                   deltaL = post_deltaL[sample],
                                                   deltaR = post_deltaR[sample],
                                                   tauL = post_tauL[sample],
                                                   tauR = post_tauR[sample])
      if (progress) {
        if(sample %% 100 == 0) {
          pb$tick(100)
        }
      }
    }
    if (method == "HPDI" ) {
      y_CI <- y_post %>%
      t(.) %>%
      coda::as.mcmc(.) %>%
      coda::HPDinterval(obj = ., prob = prob)
    }
    if (method == "ET") {
      tail <- (1 - prob) / 2
      y_CI <- y_post %>%
        apply(., 1, FUN = quantile, probs = c(0 + tail, 1 - tail)) %>%
        t(.)
    }

    result <- data.frame(transectDist = distance,
                         p_mean = y_mean,
                         p_median = y_median,
                         lower = y_CI[,1],
                         upper = y_CI[,2],
                         stringsAsFactors = F, row.names = NULL)
    names(result) <- c("transectDist", "p_mean", "p_median", low.name, up.name)

  } else {
    result <- data.frame(transectDist = distance,
                         p_mean = y_mean,
                         p_median = y_median,
                         stringsAsFactors = F)
  }

  result
})
