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

test_fit_cline <- function(stan_data, init_list, chains, model) {
  index <- which(names(stanmodels) == model)
  rstan::sampling(object = stanmodels[[index]], data = stan_data, chains = chains, init = init_list)
}

