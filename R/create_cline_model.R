#' Create the stan model code for a cline analysis
#'
#' Generate Stan model code for the specified model(s). Incorporates the priors
#' provided in \code{prior_file} into the Stan model code.
#'
#' This function is generally meant to be run internally, within the
#' \code{\link{fit_cline}} function. In that case, it is called with only one
#' option each for \code{type}, \code{tails}, and \code{direction} such that it
#' will return a list with one element containing model code for one model.
#'
#' However, this function is written flexibly such that multiple options can be
#' supplied to \code{type}, \code{tails}, and \code{direction}. In that case,
#' the returned value is a named list that contains all possible models that can
#' be made with the supplied options, see examples.
#'
#' Models are stitched together from the Stan model code provided in
#' \code{model_pieces.R}.
#'
#' @importFrom magrittr "%>%"
#'
#' @param prior_file The path to the \code{.yaml} file which contains the
#'   specifications of the priors
#' @param type The type of model to generate. Either "bi", for a binomial model
#'   of allele frequencies, or "multi" for a multinomial model of genotype
#'   frequencies.
#' @param tails Which type of tails for the model: "none", "left", "right", "mirror", or
#'   "ind"?
#' @param direction Should the model be for a cline which is increasing in
#'   frequency ("inc"), or decreasing in frequency ("dec")?
#'
#' @return A named list containing the Stan model code for the model(s) desired.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate a single multinomial model for an increasing cline with mirrored tails
#' models <- create_cline_model("prior_file.yaml", type = "multi", tails = "mirror", direction = "inc")
#' length(models)
#'
#' # Generate all 5 possible tail models for a decreasing cline with binomial likelihood
#' models <- create_cline_model("prior_file.yaml", type = "bi", direction = "dec")
#' length(models)
#' }
#'

create_cline_model <- function(prior_file,
                               type = c("bi", "multi"),
                               tails = c("none", "left", "right", "mirror", "ind"),
                               direction = c("inc", "dec")) {
  type <- match.arg(type, several.ok = T)
  tails <- match.arg(tails, several.ok = T)
  direction <- match.arg(direction, several.ok = T)

  # Set up the list that will hold all the model results
  possibilites <- expand.grid(type, tails, direction)
  model_names <- paste(possibilites[,1], possibilites[,2], possibilites[,3], sep = "_")

  result_models <- as.list(rep("NULL", times = length(model_names)))
  names(result_models) <- model_names

  # Load in the priors from the prior file
  priors <- parse_prior_file(prior_file)

  # Turn the yaml file parameters into proper lines of stan code to be pasted in to the models
  # A possible hack to get around the undefined global variable issue with R CMD Check:
  p.center <- p.deltaL <- p.deltaM <-p.deltaR <- p.f <- NULL
  p.pmax <-p.pmin <- p.tauL <- NULL
  p.tauM <- p.tauR <-p.width <- NULL
  # Seems to have worked.

  i <- 1
  for (i in 1:length(names(priors))) {
    param <- names(priors)[i]
    assign(x = paste0("p.", param), value = paste(param, " ~ ", priors[[i]], ";\n", sep = ""))
    i <- i + 1
  }

  # Combine things a bit, to make the later steps slightly less painful.
  priors.all <- paste(p.center, p.width, p.pmin, p.pmax, sep = "")
  priors.left <- paste(p.deltaL, p.tauL)
  priors.right <- paste(p.deltaR, p.tauR)
  priors.mirror <- paste(p.deltaM, p.tauM)
  priors.f <- p.f

  # Structure of Stan file: data, parameters, transformed parameters, model, generated quantities

  for (index in 1:length(result_models)) {
    if (names(result_models)[index] == "bi_none_inc") {
      result_models[index] <- paste(bi_none_inc_before_priors,
                                    priors.all,  bi_none_inc_after_priors, sep = "")
    }
    if (names(result_models)[index] == "bi_left_inc") {
      result_models[index] <- paste(bi_left_inc_before_priors,
                                    priors.all, priors.left,
                                    bi_left_inc_after_priors, sep = "")
    }
    if (names(result_models)[index] == "bi_right_inc") {
      result_models[index] <- paste(bi_right_inc_before_priors,
                                    priors.all, priors.right,
                                    bi_right_inc_after_priors, sep = "")
    }
    if (names(result_models)[index] == "bi_mirror_inc") {
      result_models[index] <- paste(bi_mirror_inc_before_priors,
                                    priors.all, priors.mirror,
                                    bi_mirror_inc_after_priors, sep = "")
    }
    if (names(result_models)[index] == "bi_ind_inc") {
      result_models[index] <- paste(bi_ind_inc_before_priors,
                                    priors.all, priors.left,
                                    priors.right, bi_ind_inc_after_priors, sep = "")
    }
    if (names(result_models)[index] == "multi_none_inc") {
      result_models[index] <- paste(multi_none_inc_before_priors,
                                    priors.all, priors.f, multi_none_inc_after_priors, sep = "")
    }
    if (names(result_models)[index] == "multi_left_inc") {
      result_models[index] <- "NEED TO MAKE"
    }
    if (names(result_models)[index] == "multi_right_inc") {
      result_models[index] <- "NEED TO MAKE"
    }
    if (names(result_models)[index] == "multi_mirror_inc") {
      result_models[index] <- "NEED TO MAKE"
    }
    if (names(result_models)[index] == "multi_ind_inc") {
      result_models[index] <- "NEED TO MAKE"
    }
  }

  result_models
}
