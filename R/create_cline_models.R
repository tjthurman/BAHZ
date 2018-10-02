#' Create the cline models
#'
#' DESCRIPTION TO BE ADDED
#'
#' DETAILS TO BE ADDED.
#'
#' @importFrom tidyr "unite"
#' @importFrom magrittr "%>%"
#'
#' @param prior_file The path to the \code{.yaml} file which contains the
#'   specifications of the priors
#'
#' @return TO BE WRITTEN
#'
#' @export
#'
#' @examples
#' #TO BE ADDED
#'

create_cline_models <- function(prior_file) {
  # Set up the list that will hold all the model results
  model_names <- expand.grid(c("bi", "multi"),
                             c("none", "left", "right", "mirror", "ind"),
                             c("inc", "dec")) %>%
    tidyr::unite(c(Var1, Var2, Var2)) %>%
    as.matrix(.)
  result_models <- as.list(rep("NULL", times = length(model_names)))
  names(result_models) <- model_names

  # Load in the priors from the prior file
  priors <- yaml.load_file(prior_file, as.named.list = T)
  # ADD a check to make sure it is 11 long, and all the right names are there

  assertthat::assert_that(length(priors) == 11, msg = "Improper number of priors, there should be 11!\nDouble-check your prior file")
  name.check <- names(priors) %in% c("center", "width", "pmin", "pmax",
                                     "deltaL", "deltaR", "deltaM",
                                     "tauL", "tauR", "tauM", "f")

  if (sum(name.check) != 11) {
    offenders <- as.vector(names(priors)[which(name.check == F)])
    stop(paste("\n", toString(offenders), "\nis/are not valid parameter names.\nDouble-check your prior file", sep = ""))
  }

  # Turn the yaml file parameters into proper lines of stan code to be pasted in to the models
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
  result_models$bi_none_inc <- paste(bi_none_inc_before_priors, "}\n",
                                    priors.all, bi_none_inc_after_priors, sep = "")

  result_models$bi_left_inc <- paste(bi_left_inc_before_priors, "}\n",
                                     priors.all, bi_left_inc_after_priors, sep = "")

  return(result_models)
}
create_cline_models("prior_config_template.yaml")
