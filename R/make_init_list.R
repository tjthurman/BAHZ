#
#' Make a list of values for initializing the Stan model.
#'
#' Generates a list, of length \code{chains}, where each element of the list is
#' a named list giving the initialization values for each model parameter.
#' Initialization values are randomly drawn from the prior distribution, and are
#' drawn independently for each chain.
#'
#' Currently supported distributions are: normal, uniform, exponential.
#'
#'
#' This function calls \code{\link{init_single_chain}} for each chain.
#' \code{\link{init_single_chain}} does all the work of making the expression
#' which will yield the proper starting values for a chain of this model.
#'
#' @export
#'
#' @param prior_file Path to the yaml file containing the prior specifications.
#'
#' @param tails Which type of tails for the model: "none", "left", "right", "mirror", or "ind."
#'
#' @param chains The number of chains to be used in Stan. Integer.
#'
#' @return A list of lists containing initialization values for the Stan model.
#'
#' @seealso \code{\link{init_single_chain}}
#'

make_init_list <- function(prior_file, tails, chains) {
  assertthat::assert_that(is.integer(chains) == T, msg = "chains must be an integer")
  assertthat::assert_that(length((tails)) == 1,
                          msg = "You must pick a single option for tails.")
  assertthat::assert_that(length((chains)) == 1,
                          msg = "You must pick a single value for chains.")

  single.chain <- init_single_chain(prior_file, tails = tails)
  init.list <- list()
  for(i in 1:chains) {
    init.list[[i]] <- eval(single.chain)
  }
  init.list
}



