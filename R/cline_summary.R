#' Summarize cline model results
#'
#' DESCRIPTION TO BE ADDED
#'
#' DETAILS TO BE ADDED. LINK TO RSTAN SUMMARY METHOD, CODA HPDI, MAYBE SOME
#' PAPERS ON WHAT THE THINGS ARE.
#'
#'
#' @param stanfit A stanfit object holding your model results.
#'
#' @param prob The probability interval to return. Default is .95. Numeric,
#'   between 0 and 1.
#'
#' @param method The method for calculating credible intervals. Either "ET" for
#'   equal-tail probability intervals, or "HPDI" for highest posterior density
#'   intervals.
#'
#' @param all By default, function returns summaries only for the main cline
#'   parameters. Set all to T
#'
#' @return A summary data frame with the columns:
#'     \itemize{
#'     \item param:
#'     \item mean:
#'     \item se_mean:
#'     \item sd:
#'     \item lower:
#'     \item upper:
#'     \item n_eff:
#'     \item Rhat:
#' }
#'
#' @export
#'
#' @examples
#' # TO ADD
#'

cline_summary <- function(stanfit, prob = .95, method = "ET", all = F) {
  # assert stanfit is a stanfit
  # assert prob is a single numeric between 0 and 1
  # assert method is ET or HPDI
  # assert all is T/F

  # Give a method to do HPDI vs. eqaul-tail interval
  tail <- (1-prob)/2
  low.name <- paste("low", prob, method, sep = "_")
  up.name <- paste("up", prob, method, sep = "_")

  # could add [abcdeghijklmnopqrstuvwxyz] in the reg expression below to
  # also keep the column with inbreeding values
  if (all == F) {
    keep <- grep("\\[|_", names(stanfit), invert = T, value = T)
  } else {
    keep <- names(stanfit)
  }

  res <- summary(stanfit, probs = c(0 + tail, 1 - tail), pars = keep)$summary %>%
    as.data.frame(.) %>%
    round(., digits = 2) %>%
    dplyr::mutate(n_eff = as.integer(.data$n_eff)) %>%
    cbind(keep, .)

  if (method == "HPDI") {
    hpd_cols <- as.matrix(stanfit, pars = keep) %>%
      coda::as.mcmc(.) %>%
      coda::HPDinterval(obj = ., prob = prob) %>%
      round(., digits = 2) %>%
      as.matrix(.)
    res[,5:6] <- hpd_cols
  }

  names(res)[c(1,5,6)] <- c("param", low.name, up.name)

  res
}
