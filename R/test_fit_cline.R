#' Fit a cline model to your data using Stan
#'
#' DESCRIPTION TO BE ADDED
#'
#' DETAILS TO BE ADDED.
#'
#'
#'
#' @export
#'

test_fit_cline <- function(stan_data, init_list, chains) {
  rstan::sampling(object = stanmodels$binom_free_none, data = stan_data, chains = chains, init = init_list)
}

