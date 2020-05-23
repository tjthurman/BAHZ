#' Summarize cline model results
#'
#' Summarizes the results of the cline model by providing the posterior mean and median,
#' credible intervals, and diagnostics for the model parameters. To save computation time, this function uses the
#' \code{\link[memoise]{memoise}} package to save past results and avoid recalculations. To clear the cache
#' of saved results, use clear.cache = T or run the \code{\link{clear_bahz_cache}} function.
#'
#' Uses the \code{rstan} summary method on class \code{\linkS4class{stanfit}} objects
#' for calculating posterior means and medians, SEM, standard deviations, equal tail
#' probability intervals, and diagnostics. Uses the \code{\link[coda]{HPDinterval}} function in
#' \code{coda} to calculate HPDI intervals, the default.
#'
#' @importClassesFrom rstan stanfit
#'
#' @import memoise
#'
#' @param stanfit A \code{\linkS4class{stanfit}} object holding your model results.
#'
#' @param prob The probability interval to return. Default is .95. Numeric,
#'   between 0 and 1.
#'
#' @param method The method for calculating credible intervals. Either "ET" for
#'   equal-tail probability intervals, or "HPDI" for highest posterior density
#'   intervals. Default is "HPDI".
#'
#' @param show.all By default, function returns summaries only for the main cline
#'   parameters. Set show.all to T to see summaries for matrix and vector parameters.
#'
#' @param clear.cache Clear the cache of saved results to ensure recalculation of summary?
#' TRUE or FALSE, default FALSE.
#'
#' @return A summary data frame with the columns:
#'     \itemize{
#'     \item param: The model parameter.
#'     \item mean: The mean of the posterior distribution.
#'     \item median: The median of the posterior distribution.
#'     \item se_mean: The standard error of the mean of the posterior distribution.
#'     \item sd: The standard deviation of the posterior distribution.
#'     \item lower: The lower limit of the credible interval. The column name
#'     will include the probability value and method used for calculating the
#'     CI.
#'     \item upper: The upper limit of the credible interval. The column name
#'     will include the probability value and method used for calculating the
#'     CI.
#'     \item n_eff: The effective number of samples from the posterior distribution.
#'     \item Rhat: The Gelman-Rubin convergence diagnostic. Should be 1, any
#'     other value is indicative of a problem during sampling.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cline_summary(stanfit)
#' }
#'

cline_summary <- memoise::memoise(function(stanfit, prob = .95, method = "HPDI",
                                           show.all = F,
                                           clear.cache = F) {
  # A bunch of argument checking
  assertthat::assert_that(class(stanfit)[1] == "stanfit",
                          msg = "Object to be summarized must be of class stanfit")
  assertthat::assert_that(is.numeric(prob) == T, msg = "prob must be numeric")
  assertthat::assert_that(length(prob) == 1, msg = "prob must be of length 1")
  assertthat::assert_that(prob <= 1, msg = "prob must be between 0 and 1")
  assertthat::assert_that(prob > 0, msg = "prob must be between 0 and 1")
  assertthat::assert_that((method %in% c("HPDI", "ET")) == T,
                          msg = "method must be either 'HPDI' or 'ET")
  assertthat::assert_that(is.logical(show.all) == T, msg = "show.all must be either TRUE or FALSE")
  assertthat::assert_that(is.logical(clear.cache) == T, msg = "clear.cache must be either TRUE or FALSE")

  if (clear.cache) {
    memoise::forget(bahz::cline_summary)
  }

  # make col names for CI columns
  tail <- (1 - prob) / 2
  low.name <- paste("low", prob, method, sep = "_")
  up.name <- paste("up", prob, method, sep = "_")

  # could add [abcdeghijklmnopqrstuvwxyz] in the reg expression below to
  # also keep the column with inbreeding values
  if (show.all == F) {
    keep <- grep("\\[|__|dev",
                 names(stanfit),
                 invert = T,
                 value = T)
  } else {
    keep <- names(stanfit)
  }

  # Initial summary, with ET probability intervals in pos
  res <-
    rstan::summary(
      stanfit,
      probs = c(0 + tail, 0.5, 1 - tail),
      pars = keep,
      use_cache = F)$summary %>%
    as.data.frame(.) %>%
    round(., digits = 2) %>%
    dplyr::mutate(n_eff = as.integer(.data$n_eff)) %>%
    dplyr::select(mean, median = `50%`, dplyr::everything()) %>%
    cbind(keep, .)

  # If doing HPDI, make a matrix of those columns and replace the ET
  # columns from the current results
  if (method == "HPDI") {
    hpd_cols <- as.matrix(stanfit, pars = keep) %>%
      coda::as.mcmc(.) %>%
      coda::HPDinterval(obj = ., prob = prob) %>%
      round(., digits = 2) %>%
      as.matrix(.)
    res[, 6:7] <- hpd_cols
  }

  # rename columns
  names(res)[c(1, 6, 7)] <- c("param", low.name, up.name)

  res
})
