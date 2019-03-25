#' Make the expression of init values for a single chain in Stan
#'
#' This function it a bit complicated. It reads in the \code{prior_file}, and
#' parses both (1) the prior distribution for a parameter, and (2) the value(s)
#' which describe that distribution. From that information, it creates for each
#' parameter an \code{rlang} expression of the appropriate R distribution
#' function (\code{rnorm}, \code{runif}, etc.) that will generate a single random
#' starting value from the prior. Finally, it assembles the expressions for each
#' parameter into a a single list using only the parameters needed for the
#' chosen model as specified by the \code{tails} argument. That single list will
#' eventually be evaluated by the \link{prep_init_list} function to create the
#' list of initial values for Stan.
#'
#' Currently supported distributions are: normal, uniform, exponential.
#'
#' Used internally, in \link{prep_init_list}.
#'
#' @keywords internal
#'
#' @importFrom stats "rexp" "rnorm" "runif"
#'
#' @param prior_file Path to the yaml file containing the prior specifications.
#'
#' @param tails Which type of tails for the model: "none", "left", "right", "mirror", or "ind"?
#'
#' @return An \code{rlang} expression of the initialization values for a single chain
#'
#'

init_single_chain <- function(prior_file, tails = c("none", "left", "right", "mirror", "ind")) {
  # Argument check
  tails <- match.arg(tails, several.ok = F)
  # Initialize all the objects I'll be making later to avoid R CMD CHECK warnings.
  init.center <- init.width <- init.pmin <- init.pmax <- NULL
  init.f <- init.deltaL <- init.deltaR <- init.deltaM <- NULL
  init.tauL <- init.tauR <- init.tauM <- NULL

  prior_list <- parse_prior_file(prior_file)

  # Originally tried to pull the for loop below into its own function,
  # parse_distribution.
  # But, the assign commands put everything into the global environment,
  # and didn't keep them here within the init.chain function.
  # Not sure why that was yet, something with the environments that I don't yet understand.
  # But keeping it here fixed the problem,
  # the inits don't get assigned to the global environment.
  # So a little ugly, but let's just leave it as is.

  for (i in 1:length(prior_list)) { # for each parameter
    # find which distribution it is
    # Based on the distribution, extract the proper values
    # then assemble that into an expression for the parameter
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

  # Assemble all the per-parameter expressions into one list for a single chain
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
    # Shouldn't really get here, with the argument checking.
    stop("Tails argument not properly specified.\nMust be either:\n`none`, `left`, `right`, `mirror`, `ind`")
  }
  init.chain
}
