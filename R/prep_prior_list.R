#' Prepare the list of prior values to be fed into Stan
#'
#' Parses the prior specifications provided in the \code{.yaml}
#' file, including checks that they are properly specified and
#' currently supported. Then, converts the prior specifications into
#' a named list to be passed to Stan.
#'
#' @param prior_file Path to the \code{.yaml} file containing the prior specifications.
#'
#' @return A list containing the values for the priors to be passed to Stan.
#'
#' @export
#'
#' @seealso \code{\link{fit_geno_cline}}
#'
#' @examples
#' \dontrun{
#' prep_prior_list("path/to/priors.yaml")
#' }


prep_prior_list <- function(prior_file) {

  # Read in the config file, which checks for proper order.
  prior_config <- parse_prior_file(prior_file)

  # Get the disributions specified for each parameter
  specified.dist <- rep(as.character(NA), times = length(prior_config))
  for (i in 1:length(prior_config)) {
    specified.dist[i] <- stringr::str_extract(prior_config[[i]], "^[a-z]*\\(") %>%
      stringr::str_remove_all("\\(")
  }

  # Then some checks
  dist_stan_integers <- rep(as.integer(NA), times = length(prior_config))
  for (i in 1:length(prior_config)) {
    # Check if the distribution is supported
    if (check_prior_supported(names(prior_config)[i], specified.dist[i]) == F) {
      stop(paste("\nThe distribution specified for ", names(prior_config)[i], " is not currently supported",
                 "\nPick a different distribution, or submit a feature request on Github", sep = ""))
    }
    # Check if the distribution is properly specified
    if (check_prior_specification(specified.dist[i], prior_config[[i]]) == F) {
      stop(paste("\nThe specification of values for ", names(prior_config)[i], " is not correct",
                 "\nCheck that there are the proper number of comma-separated values for your distribution",
                 "\n E.g., 2 values for a normal distribution, 1 value for an exponential, etc.", sep = ""))
    }
    # And generate the STAN integer in put for prior types
    dist_stan_integers[i] <- assign_stan_dist_int(specified.dist[i])
  }

  # Turn that into a list of the distribution choices
  df <- data.frame(dist_stan_integers)
  df$stan_prior_names <- c("p_dist_center",
                           "p_dist_width",
                           "p_dist_pmin",
                           "p_dist_pmax",
                           "p_dist_deltaL",
                           "p_dist_deltaR",
                           "p_dist_deltaM",
                           "p_dist_tauL",
                           "p_dist_tauR",
                           "p_dist_tauM",
                           "p_dist_f")
  prior.type <- df %>%
    dplyr::select(.data$stan_prior_names, .data$dist_stan_integers) %>%
    tidyr::spread(.data$stan_prior_names, .data$dist_stan_integers) %>%
    as.list(.)

  # Now, need to extract the values of the distributions.
  val1 <- rep(999, length(prior_config))
  val2 <- rep(999, length(prior_config))
  for (i in 1:length(prior_config)) {
    if (specified.dist[i] %in% c("normal", "uniform", "beta")) {
      val1[i] <- as.numeric(extract_first(prior_config[[i]]))
      val2[i] <- as.numeric(extract_last(prior_config[[i]]))
    }
    if (specified.dist[i] %in% c("exponential", "poisson")) {
      val1[i] <- as.numeric(extract_only(prior_config[[i]]))
    }
  }

  # Turn that intp a list of prior specifications
  df2 <- data.frame(param = names(prior_config),
             distribution = specified.dist,
             val1 = val1,
             val2 = val2) %>%
    tidyr::gather(.data$val1, .data$val2, key = "which", value = "value") %>%
    dplyr::arrange(match(.data$param, names(prior_config)))

  df2$stan_variable_names <- c("p_center_1", "p_center_2",
                               "p_width_1", "p_width_2",
                               "p_min_1", "p_min_2",
                               "p_max_1", "p_max_2",
                               "p_deltaL_1", "p_deltaL_2",
                               "p_deltaR_1", "p_deltaR_2",
                               "p_deltaM_1", "p_deltaM_2",
                               "p_tauL_1", "p_tauL_2",
                               "p_tauR_1", "p_tauR_2",
                               "p_tauM_1", "p_tauM_2",
                               "p_f_1", "p_f_2")

  prior.value <- df2 %>%
    dplyr::select(.data$stan_variable_names, .data$value) %>%
    tidyr::spread(.data$stan_variable_names, .data$value) %>%
    as.list(.)
  # And combine the lists together
  c(prior.type, prior.value)
}

