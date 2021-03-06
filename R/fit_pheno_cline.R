#' Fit a phenotypic cline model to your data
#'
#' Use Rstan to fit Bayesian hybrid zone cline models to phenotypic data.
#' The phenotypic cline models in BAHZ assume that phenotypes are normally distributed.
#'
#' This is a wrapper function, which calls various data- and model-preparation
#' functions from \code{bahz} before passing the results to the
#' \code{\link[rstan]{sampling}} function in \code{rstan}, which does the MCMC model
#' fitting. Specifically, this function calls:
#'
#' \enumerate{
#'  \item \code{\link{prep_pheno_data}} to format the provided
#'  \code{data}.
#'  \item \code{\link{prep_prior_list}}, which creates a list containing
#'  the specifications of the prior distributions, to be passed to Stan.
#'  It creates the priors according to the configuration file \code{prior_file}.
#'  \item \code{\link{prep_init_list}}, which creates a list of values for
#'  initialiing the MCMC chains in \code{stan}. Initial values are randomly
#'  drawn from the priors specified in \code{prior_file}. This step is ignored
#'  if the user supplies an initialization list with the \code{init} argument.
#' }
#'
#' The results of these functions are then passed to
#' \code{\link[rstan]{sampling}}, along with any additional \code{stan}
#' arguments supplied by the user. In most cases, these additional \code{stan}
#' arguments should be unnecessary.  \code{bahz} uses the default \code{stan}
#' arguments for all sampler parameters expect for one: the default adapt_delta
#' for \code{bahz} is 0.95, to better avoid divergent transitions. Note that
#' supplying ANY \code{stan} arguments will override ALL \code{bahz} defaults.
#' E.g., specifying a different number of total iterations will reset the
#' adapt_delta to the \code{stan} default of 0.8.
#'
#' @param data A dataframe containing your cline data. See
#'   \code{\link{prep_pheno_data}} for possible formats.
#' @param prior_file The path to the \code{.yaml} file which contains the
#'   specifications of the priors.
#' @param pheno_variance The method used to estimate phenotypic variance.
#'   Either "constant" (the default), "pooled", or "independent". If "constant", bahz
#'   estimates a single, constant phenotypic variance across all sites
#'   (global_sigma). If "independent", bahz estimates phenotypic variance
#'   independently at each site (site_sigma). If "pooled", bahz uses a
#'   hierarchical model which shares information across sites to estimate
#'   both the average phenotypic variance (global_sigma) and site-specific
#'   variances (site_sigma). The site-specific variances are assumed to be
#'   normally distributed around the global_sigma with their own variance,
#'   which is estimated (sigma_sigma).
#' @param chains The number of MCMC chains to create. Numeric, coerced to
#'   integer. Default is 4.
#' @param init Optional, default is \code{NULL}. A user-provided list of
#'   initialization values for \code{stan}, to be used instead of the
#'   \code{bahz} default of random initialization values from the prior. See
#'   \code{\link[rstan]{stan}} for details on how to specify the init list.
#' @param ignore_data Used for doing prior predictive checks. If FALSE, the default,
#'   the model is fit normally. If TRUE, the data are ignored during model fitting.
#'   This is useful for double-checking the prior distributions or parameters, and
#'   for doing prior predictive checks of the model.
#' @param ... Arguments to be passed to \code{stan}, e.g., number of iterations,
#'   warmup period, etc. See \code{\link[rstan]{sampling}} and the
#'   \code{control} argument in \code{\link[rstan]{stan}} for information on
#'   these arguments, and the details below for information on how bahz defaults
#'   interact with stan defaults.
#'
#' @return A \code{\linkS4class{stanfit}} object containing your model results.
#'
#' @export
#'
#'
#' @seealso \code{\link{prep_pheno_data}}, \code{\link{prep_prior_list}},
#'   \code{\link{prep_init_list}}, \code{\link[rstan]{sampling}}
#'
#' @examples
#' \dontrun{
#' # Fit a phenotypic cline
#' # Uses default number of chains and stan parameters.
#' results <- fit_pheno_cline(clinedata, "prior_file.yaml")
#'
#' # Fit a phenotypic cline
#' # Use 5 chains with 4000 warmup iterations
#' # and 8000 total iterations per chain.
#' # Note that this changes the adapt_delta to
#' # the stan default of 0.8
#'
#' results2 <- fit_pheno_cline(clinedata, "prior_file.yaml",
#'                      chains = 5,
#'                      iter = 8000, warmup = 4000)
#' }
#'


fit_pheno_cline <- function(data, prior_file, pheno_variance = "constant",
                           chains = 4,
                           init = NULL, ignore_data = F, ...) {
  # Argument checking
  assertthat::assert_that(is.numeric(chains) == T, msg = "chains must be numeric")
  pheno_variance <- match.arg(pheno_variance, choices = c("constant", "pooled", "independent"), several.ok = F)
  assertthat::assert_that(is.logical(ignore_data), msg = "ignore_data must be either TRUE or FALSE")
  ch <- as.integer(chains)

  # Calling internal functions
  # Prep data
  stan_data <- prep_pheno_data(data)
  # Make list of prior values
  # this also runs a bunch of prior compatibility checks
  prior_list <- prep_prior_list(prior_file)
  # Make list of initial values
  if (is.null(init)) {
    init_list <- prep_init_list(prior_file, tails = "none", chains = ch, type = "pheno")
  } else {
    init_list <- init
  }

  # Process ignore data
  if (ignore_data) {# if ignoring the data
    ignoreData <- list(ignoreData = as.integer(1))
  } else {# otherwise, use the data
    ignoreData <- list(ignoreData = as.integer(0))
  }

  # Find location of the model in the stanmodels object that matches
  # the desired model provide by the user
  model_index <- which(names(stanmodels) == paste("pheno", pheno_variance, sep = "_"))

  # Pass everything to stan
  if (length(eval(substitute(alist(...)))) > 0) {# if user supplies extra parameters to go to Stan
    clinefit <- rstan::sampling(object = stanmodels[[model_index]], data = c(stan_data, prior_list, ignoreData),
                                chains = ch, init = init_list, ...)
  } else {# otherwise, use the bahz defaults
    clinefit <- rstan::sampling(object = stanmodels[[model_index]], data = c(stan_data, prior_list, ignoreData),
                                chains = ch, init = init_list, control = list(adapt_delta = 0.95))
  }

  clinefit
}
