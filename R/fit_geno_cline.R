#' Fit a genetic cline model to your data
#'
#' Use Rstan to fit Bayesian hybrid zone cline models to genetic data.
#'
#' This is a wrapper function, which calls various data- and model-preparation
#' functions from \code{bahz} before passing the results to the
#' \code{\link[rstan]{sampling}} function in \code{rstan}, which does the MCMC model
#' fitting. Specifically, this function calls:
#'
#' \enumerate{
#'  \item \code{\link{prep_geno_data}} to format the provided
#'  \code{data}.
#'  \item \code{\link{prep_prior_list}}, which creates a list containing
#'  the specifications of the prior distributions, to be passed to Stan.
#'  It creates the priors according to the configuration file \code{prior_file}.
#'  \item \code{\link{prep_init_list}}, which creates a list of values for
#'  initialiing the MCMC chains in \code{stan}. Initial values are randomly
#'  drawn from the priors specified in \code{prior_file}.
#' }
#'
#' The results of these functions are then passed to \code{\link[rstan]{sampling}}, along with
#' any of the additional \code{stan} arguments supplied by the user.
#'
#' @param data A dataframe containing your cline data. See
#'   \code{\link{prep_geno_data}} for possible formats.
#' @param prior_file The path to the \code{.yaml} file which contains the
#'   specifications of the priors.
#' @param type The type of model to fit Either "bi", for a binomial model
#'   of allele frequencies, or "multi" for a multinomial model of genotype
#'   frequencies.
#' @param tails Which type of tails for the model: "none", "left", "right", "mirror", or
#'   "ind"?
#' @param chains The number of MCMC chains to create. Numeric, coerced to
#'   integer. Default is 3.
#' @param ... Arguments to be passed to \code{stan}, e.g., number of iterations, warmup
#'   period, etc. See \code{\link[rstan]{sampling}} for more information.
#'
#' @return A \code{\linkS4class{stanfit}} object containing your model results.
#'
#' @export
#'
#'
#' @seealso \code{\link{prep_geno_data}}, \code{\link{prep_prior_list}},
#'   \code{\link{prep_init_list}}, \code{\link[rstan]{sampling}}
#'
#' @examples
#' \dontrun{
#' # Fit a multinomial cline
#' # with mirrored introgression tails.
#' # Uses default number of chains and stan parameters.
#' results <- fit_geno_cline(clinedata, "prior_file.yaml",
#'                      type = "multi", tails = "mirror")
#'
#' # Fit a binomial cline
#' # with no introgression tails.
#' # Use 4 chains with 4000 warmup iterations
#' # and 8000 total iterations per chain.
#'
#' results2 <- fit_geno_cline(clinedata, "prior_file.yaml",
#'                      type = "bi", tails = "none",
#'                      chains = 4,
#'                      iter = 8000, warmup = 4000)
#' }
#'
#'
#'
#'

fit_geno_cline <- function(data, prior_file,
                      type = c("bi", "multi"),
                      tails = c("none", "left", "right", "mirror", "ind"),
                      chains = 3, ...) {
  # Argument checking
  type <- match.arg(type, several.ok = F)
  tails <- match.arg(tails, several.ok = F)
  assertthat::assert_that(is.numeric(chains) == T, msg = "chains must be numeric")
  ch <- as.integer(chains)

  # Calling internal functions
  # Prep data
  stan_data <- prep_geno_data(data, type = type)
  # Make list of prior values
  # this also runs a bunch of prior compatibility checks
  prior_list <- prep_prior_list(prior_file)
  # Make list of initial values
  init_list <- prep_init_list(prior_file, tails = tails, chains = ch)

  # Find location of the model in the stanmodels object that matches
  # the desired model provide by the user
  model.index <- which(names(stanmodels) == paste(type, "free", tails, sep = "_"))

  # Pass everything to stan
  clinefit <- rstan::sampling(object = stanmodels[[model.index]], data = c(stan_data, prior_list),
                              chains = ch, init = init_list, ...)

  clinefit
}
