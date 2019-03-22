#
#' Make a list of vlaues for initializing the stan model.
#'
#' This function calls \link{init_single_chain}, which does all the work of
#' making the expression which will yeild the proper starting values for a chain
#' of this model. It then evaluates that expression for each desired chain in
#' the model and adds it to a list.
#'
#' Currently supported distributions are: normal, uniform, exponential.
#'
#' @export
#'
#' @param prior_file Path to the yaml file containing the prior specifications.
#'
#' @param tails Which type of tails for the model: "none", "left", "right", "mirror", or "ind."
#'
#' @param chains The number of chains to be used in Stan.
#'
#' @return A list of lists containing initialization values for the Stan model.
#'
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
  # INSERT SOME CHECKs HERE OF THE VALUES.
  # Specifically, make sure width is >0.
  # Or maybe just add abs() for width in the init_single_chain fxn?
  init.list
}



