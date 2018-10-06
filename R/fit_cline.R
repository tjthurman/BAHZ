#' Fit a cline model to your data using Stan
#'
#' DESCRIPTION TO BE ADDED
#'
#' DETAILS TO BE ADDED.
#'
#'
#' @param data TO BE WRITTEN
#' @param prior_file TO BE WRITTEN
#' @param type TO BE WRITTEN
#' @param tails TO BE WRITTEN
#' @param direction TO BE WRITTEN
#' @param chains The number of MCMC chains to create. Integer. Default is 3.
#' @param ... Arguments to be passed to stan. LINK TO IT.
#'
#' @return A stanfit object containing your model results.
#'
#' @export
#'
#' @examples
#' #TO BE ADDED
#'

fit_cline <- function(data, prior_file, type, tails, direction, chains = 3, ...) {
  for (singular.args in alist(type, tails, direction, chains)) {
    assertthat::assert_that(length(eval(singular.args)) == 1,
                            msg = paste("You must pick a single value for", singular.args, sep = " "))
  }
  assertthat::assert_that(is.numeric(chains) == T, msg = "chains must be numeric")
  ch <- as.integer(chains)


  stan_data <- load_cline_data(data, type = type)
  stan_model <- create_cline_model(prior_file, type = type, tails = tails, direction = direction)[[1]]
  init_list <- make_init_list(prior_file, tails = tails, chains = ch)

  clinefit <- rstan::stan(model_code = stan_model, data = stan_data, chains = ch, init = init_list, ...)

  clinefit
}
