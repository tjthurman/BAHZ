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

test_fit_cline2 <- function(stan_data, init_list, chains) {
  rstan::sampling(object = stanmodels$minimal, data = stan_data, chains = chains, init = init_list)
}

