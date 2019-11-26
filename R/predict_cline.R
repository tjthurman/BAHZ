#' Predict allele frequencies or mean phenotypes for a cline
#'
#' Uses the cline parameters from your stanfit object (the mean of the posterior
#' distribtution of each parameter) to predict the allele frequency (for genetic
#' clines) or the mean phenotype value (for phenotypic clines) at the
#' distance(s) along the transect specified by the user. To save computation time, this function use the
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
#'  \item p: the predict allele frequency or mean phenotype at each distance,
#'  according to the cline equation. Numeric.
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

predict_cline <- memoise::memoise(function(stanfit, distance,
                                           confidence = F, prob = 0.95,
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

  #Pass those to general cline equation to get the
  # best fit cline
  y <- general_cline_eqn(transectDist = distance,
                         decrease = decreasing,
                         center = center,
                         width = width,
                         pmin = pmin,
                         pmax = pmax,
                         deltaL = deltaL,
                         tauL = tauL,
                         deltaR = deltaR,
                         tauR = tauR)



  # If we want to add confidence intervals.
  # Could possible speed this up a bit by writing around the loop somehow.
  # Not sure.
  if (confidence == T) {
    # make col names for CI columns
    # tail <- (1 - prob) / 2
    low.name <- paste("low", prob, "HPDI", sep = "_")
    up.name <- paste("up", prob, "HPDI", sep = "_")

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
        post_deltaR <- NULL
        post_tauR <- NULL
      }
    } else if ("deltaR" %in% colnames(posterior)) { # if right only
      post_deltaR <-  posterior[ , which(colnames(posterior) == "deltaR")]
      post_tauR <- posterior[ , which(colnames(posterior) == "tauR")]
      post_deltaL <- NULL
      post_tauL <- NULL
    } else { # no tails
      post_deltaL <- NULL
      post_deltaR <- NULL
      post_tauL <- NULL
      post_tauR <- NULL
    }
    if (progress) {
      pb <- progress::progress_bar$new(total = dim(posterior)[1],
                             format = " Posterior sample [:bar] :current/:total",
                             clear = F, width= 60)
    }

    pb$tick(0)
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
    y_CI <- y_post %>%
      t(.) %>%
      coda::as.mcmc(.) %>%
      coda::HPDinterval(obj = ., prob = prob)

    result <- data.frame(transectDist = distance,
                         p = y,
                         lower = y_CI[,1],
                         upper = y_CI[,2],
                         stringsAsFactors = F, row.names = NULL)
    names(result) <- c("transectDist", "p", low.name, up.name)

  } else {
    result <- data.frame(transectDist = distance,
               p = y, stringsAsFactors = F)
  }

  result
})
