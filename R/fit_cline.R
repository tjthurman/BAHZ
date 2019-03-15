#' Fit a cline model to your data
#'
#' DESCRIPTION TO BE ADDED
#'
#' This is a wrapper function, which calls various data- and model-preparation
#' functions from \code{bahz} before passing the results to the
#' \code{\link[rstan]{stan}} function in \code{rstan}, which does the MCMC model
#' fitting. Specifically, this function calls:
#'
#' \enumerate{
#'   \item \code{\link{load_cline_data}} to format the provided
#'   \code{data}.
#'  \item \code{\link{create_cline_model}} to create the Stan model code for the
#'  model, as specified by \code{type}, \code{tails}, \code{direction}, and the
#'  priors in \code{prior_file}.
#'  \item \code{\link{make_init_list}}, which creates a list of values for
#'  initialiing the MCMC chains in \code{stan}. Initial values are randomly
#'  drawn from the priors specified in \code{prior_file}.
#' }
#'
#' The results of these functions are then passed to \code{stan}, along with
#' any of the additional \code{stan} arguments supplied by the user.
#'
#' @param data A dataframe containing your cline data. See
#'   \code{\link{load_cline_data}} for possible formats.
#' @param prior_file The path to the \code{.yaml} file which contains the
#'   specifications of the priors
#' @param type The type of model to generate. Either "bi", for a binomial model
#'   of allele frequencies, or "multi" for a multinomial model of genotype
#'   frequencies.
#' @param tails Which type of tails for the model: "none", "left", "right", "mirror", or
#'   "ind"?
#' @param direction Should the model be for a cline which is increasing in
#'   frequency ("inc"), or decreasing in frequency ("dec")?
#' @param chains The number of MCMC chains to create. Numeric, coerced to
#'   integer. Default is 3.
#' @param ... Arguments to be passed to stan, e.g., number of iterations, warmup
#'   period. See \code{\link[rstan]{stan}} for more information.
#'
#' @return A \code{\linkS4class{stanfit}} object containing your model results.
#'
#' @export
#'
#' @seealso \code{\link{load_cline_data}}, \code{\link{create_cline_model}},
#'   \code{\link{make_init_list}}, \code{\link[rstan]{stan}}
#'
#' @examples
#' \dontrun{
#' # Fit a decreasing, multinomial cline
#' # with mirrored introgression tails.
#' # Uses default number of chains, and stan parameters.
#' results <- fit_cline(clinedata, "prior_file.yaml",
#'                      type = "multi", tails = "mirror",
#'                      direction = "dec")
#'
#' # Fit a decreasing, binomial cline
#' #  with no introgression tails.
#' # Use 4 chains with 4000 warmup iterations
#' # and 8000 total iterations per chain.
#'
#' results2 <- fit_cline(clinedata, "prior_file.yaml",
#'                      type = "bi", tails = "none",
#'                      direction = "dec", chains = 4,
#'                      iter = 8000, warmup = 4000)
#' }
#'
#'
#'
#'

fit_cline <- function(data, prior_file,
                      type = c("bi", "multi"),
                      tails = c("none", "left", "right", "mirror", "ind"),
                      direction = c("inc", "dec"),
                      chains = 3, ...) {
  type <- match.arg(type, several.ok = F)
  tails <- match.arg(tails, several.ok = F)
  direction <- match.arg(direction, several.ok = F)

  assertthat::assert_that(is.numeric(chains) == T, msg = "chains must be numeric")
  ch <- as.integer(chains)


  stan_data <- load_cline_data(data, type = type)
  stan_model <- create_cline_model(prior_file, type = type, tails = tails, direction = direction)[[1]]
  init_list <- make_init_list(prior_file, tails = tails, chains = ch)

  clinefit <- rstan::stan(model_code = stan_model, data = stan_data, chains = ch, init = init_list, ...)

  clinefit
}
