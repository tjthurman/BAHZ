#
#' Make the expression of init values for a single chain
#'
#' This function it a bit complicated. It reads in the prior_file, and parses
#' both (1) the prior distribution for a parameter, and (2) the value(s) which
#' describe that distribution. It saves that info as an expression of the
#' appropriate R function to generate a single random value from that
#' distribution, for each parameter. Then, it assembles those parameters into a
#' list, only using the parameters needed for a certain model, specified by tails.
#'
#' Currently supported distributions are: normal, uniform, exponential.
#'
#' Used internally, in \link{make_init_list}.
#'
#' @keywords internal
#'
#' @importFrom stats "rexp" "rnorm" "runif"
#'
#' @param prior_file Path to the yaml file containing the prior specifications.
#'
#' @param tails Which type of tails for the model: "none", "left", "right", "mirror", or "ind."
#'
#' @return An rlang expression of the initialization values for a single chain
#'
#'
#'

init_single_chain <- function(prior_file, tails = c("none", "left", "right", "mirror", "ind")) {
  tails <- match.arg(tails, several.ok = F)

  init.center <- init.width <- init.pmin <- init.pmax <- NULL
  init.f <- init.deltaL <- init.deltaR <- init.deltaM <- NULL
  init.tauL <- init.tauR <- init.tauM <- NULL

  prior_list <- parse_prior_file(prior_file)

  # Originally tried to pull the for loop below into its own function,
  # parse_distribution.
  # But, the assign commands put everything into the global environment,
  # and didn't keep them here within the init.chain function.
  # Not sure why that was yet, something with the environments that I don't yet understand.
  # But I keeping it here fixed the problem,
  # the inits don't get assigned to the global environment.

  for (i in 1:length(prior_list)) {
    if (stringr::str_count(prior_list[[i]], "^normal\\(") == 1) {
      mean <- extract_first(prior_list[[i]])
      sd <- extract_last(prior_list[[i]])
      init.name <- paste("init", names(prior_list)[i], sep = ".")
      assign(x = init.name, value = rlang::expr(rnorm(n = 1, mean = !!mean, sd = !!sd)), inherits = T)
    } else if (stringr::str_count(prior_list[[i]], "^uniform\\(") == 1) {
      low <- extract_first(prior_list[[i]])
      up <- extract_last(prior_list[[i]])
      init.name <- paste("init", names(prior_list)[i], sep = ".")
      assign(x = init.name, value = rlang::expr(runif(n = 1, min = !!low, max = !!up)), inherits = T)
    } else if (stringr::str_count(prior_list[[i]], "^exponential\\(") == 1) {
      rate <- extract_only(prior_list[[i]])
      init.name <- paste("init", names(prior_list)[i], sep = ".")
      assign(x = init.name, value = rlang::expr(rexp(n = 1, rate = !!rate)), inherits = T)
    } else {
      message <- paste0("The distribution you've selected for parameter\n",
                        names(prior_list[i]), "\n",
                        "Isn't supported yet. Submit a feature request on github!")
      stop(message)
    }
  }

  if (tails == "none") {
    init.chain <- rlang::expr(list(center = !!init.center,
                                   width = !!init.width,
                                   pmin = !!init.pmin,
                                   pmax = !!init.pmax))
  } else if (tails == "left") {
    init.chain <- rlang::expr(list(center = !!init.center,
                                   width = !!init.width,
                                   pmin = !!init.pmin,
                                   pmax = !!init.pmax,
                                   deltaL = !!init.deltaL,
                                   tauL = !!init.tauL))
  } else if (tails == "right") {
    init.chain <- rlang::expr(list(center = !!init.center,
                                   width = !!init.width,
                                   pmin = !!init.pmin,
                                   pmax = !!init.pmax,
                                   deltaR = !!init.deltaR,
                                   tauR = !!init.tauR))
  } else if (tails == "mirror") {
    init.chain <- rlang::expr(list(center = !!init.center,
                                   width = !!init.width,
                                   pmin = !!init.pmin,
                                   pmax = !!init.pmax,
                                   deltaM = !!init.deltaM,
                                   tauM = !!init.tauM))
  } else if (tails == "ind") {
    init.chain <- rlang::expr(list(center = !!init.center,
                                   width = !!init.width,
                                   pmin = !!init.pmin,
                                   pmax = !!init.pmax,
                                   deltaL = !!init.deltaL,
                                   tauL = !!init.tauL,
                                   deltaR = !!init.deltaR,
                                   tauR = !!init.tauR))
  } else {
    stop("Tails argument not properly specified.\nMust be either:\n`none`, `left`, `right`, `mirror`, `ind`")
  }
  return(init.chain)
}
