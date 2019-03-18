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

test_fit_cline <- function(data, init_list) {
  rstan::sampling(stanmodels$minimal, data = stan_data, chains = 3, init = init_list)
}
