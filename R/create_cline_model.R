#' Create the model code for a cline analysis
#'
#' DESCRIPTION TO BE ADDED
#'
#' DETAILS TO BE ADDED.
#'
#' @importFrom magrittr "%>%"
#'
#' @param prior_file The path to the \code{.yaml} file which contains the
#'   specifications of the priors
#' @param type TO BE WRITTEN
#' @param tails TO BE WRITTEN
#' @param direction TO BE WRITTEN
#'
#' @return TO BE WRITTEN
#'
#' @export
#'
#' @examples
#' #TO BE ADDED
#'

create_cline_model <- function(prior_file,
                               type = c("bi", "multi"),
                               tails = c("none", "left", "right", "mirror", "ind"),
                               direction = c("inc", "dec")) {
  # Set up the list that will hold all the model results
  possibilites <- expand.grid(type, tails, direction)
  model_names <- paste(possibilites[,1], possibilites[,2], possibilites[,3], sep = "_")

  result_models <- as.list(rep("NULL", times = length(model_names)))
  names(result_models) <- model_names

  # Load in the priors from the prior file
  priors <- yaml::yaml.load_file(prior_file, as.named.list = T)
  # ADD a check to make sure it is 11 long, and all the right names are there

  assertthat::assert_that(length(priors) == 11, msg = "Incorrect number of priors, there should be 11!\nDouble-check your prior file")
  name.check <- names(priors) %in% c("center", "width", "pmin", "pmax",
                                     "deltaL", "deltaR", "deltaM",
                                     "tauL", "tauR", "tauM", "f")

  if (sum(name.check) != 11) {
    offenders <- as.vector(names(priors)[which(name.check == F)])
    stop(paste("\n", toString(offenders), "\nis/are not valid parameter names.\nDouble-check your prior file", sep = ""))
  }

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
