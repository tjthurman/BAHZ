#' Make the list of prior value to be fed into Stan
#'
#' DESCRIPTION TO BE ADDED
#'
#' DETAILS TO BE ADDED.
#'
#'
#' @param prior_file TO BE WRITTEN
#'
#' @return A list containing the values for the priors to be passed to Stan
#'
#' @export
#'
#' @examples
#' #TO BE ADDED


make_prior_list <- function(prior_file) {

  # out.prior.list <- list(p_m_center = 0,
  #                        p_sd_center = 0,
  #                        p_m_width = 0,
  #                        p_sd_width = 0,
  #                        p_scale_width = 0,
  #                        p_l_min = 0,
  #                        p_u_min = 0.2,
  #                        p_l_max = 0.8,
  #                        p_u_max = 1)


  # prior_file <- "prior_config_template.yaml"
  # # Read in the config file, which checks for proper order.
  # prior_config <- parse_prior_file(prior_file)
  #
  #
  # # Get the disributions specified for each parameter
  # specified.dist <- rep(as.character(NA), times = length(prior_config))
  # for (i in 1:length(prior_config)) {
  #   specified.dist[[i]] <- stringr::str_extract(prior_config[[i]], "^[a-z]*\\(") %>%
  #     stringr::str_remove_all("\\(")
  # }
  #
  # # Then some checks
  # for (i in 1:length(prior_config)) {
  #   # Check if the distribution is supported
  #   if (check_prior_supported(names(prior_config)[i], specified.dist[i]) == F) {
  #     stop(paste("\nThe distribution specified for", names(prior_config)[i], "is not currently supported",
  #                "\nPick a different distribution, or submit a feature request on Github", sep = ""))
  #   }
  #   # Check if the distribution is properyl specified
  #   if (check_prior_specification(specified.dist[i], prior_config[[i]]) == F) {
  #     stop(paste("\nThe specification of values for", names(prior_config)[i], "is not correct",
  #                "\nCheck that there are the proper number of comma-separated values for your distribution",
  #                "\n E.g., 2 values for a normal distribution, 1 value for an exponential, etc.", sep = ""))
  #   }
  # }
  #
  # # Now, need to extract the values.
  # val1 <- rep(999, length(prior_config))
  # val2 <- rep(999, length(prior_config))
  # for (i in 1:length(prior_config)) {
  #   print(specified.dist[i])
  #   if (specified.dist[i] %in% c("normal", "uniform", "beta")) {
  #     val1[i] <- as.numeric(extract_first(prior_config[[i]]))
  #     val2[i] <- as.numeric(extract_last(prior_config[[i]]))
  #   }
  #   if (specified.dist[i] %in% c("exponential", "poisson")) {
  #     val1[i] <- as.numeric(extract_only(prior_config[[i]]))
  #   }
  # }
  # STILL NEED TO:
  # GET THE VALUE OF WHICH DISTRIBUTION TO USE (0/1/2 etc)
  # Turn all these results into a list somehow.

 # z <-  data.frame(name = c("a", "b"), value = c(1,2)) %>%
 #      spread(name, value, convert = F) %>%
 #    as.list(.)

  # out.prior.list
}

