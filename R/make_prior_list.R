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

#
#   prior_file <- "prior_config_template.yaml"
#   prior_config <- parse_prior_file(prior_file)
#
#   out.prior.list <- list(p_m_center = 0,
#                          p_sd_center = 0,
#                          p_m_width = 0,
#                          p_sd_width = 0,
#                          p_scale_width = 0,
#                          p_l_min = 0,
#                          p_u_min = 0.2,
#                          p_l_max = 0.8,
#                          p_u_max = 1)
#   for (i in 1:length(prior_config)) {
 #    param <- names(prior_config)[i]
 #    print(param)
 #    distribution <- stringr::str_extract(prior_config[[i]], "^[a-z]*\\(") %>%
 #      stringr::str_remove_all("\\(")
 #    # check compatibility
 #    # a little function that takes the parameter and the distribution as input and says if they're OK
 #    # to go in the utilities.
 #
 #    # maybe put below in a function as well?
 #    # could have a "supported_priors.csv" in the package with the skeleton of the priors
 #    system.file("extdata", "prior_config_template.yaml", package = "bahz", mustWork = TRUE),
 #    if (distribution %in% c("normal", "uniform")) {
 #      # some sort of check that there actually are 2 value? Or just rely on the checks in the extract functions?
 #      val1 <- extract_first(prior_config[[i]])
 #      val2 <- extract_last(prior_config[[i]])
 #      print(distribution)
 #      print(paste(val1, val2))
 #    } else if (distribution %in% c("exponential")) {
 #      val1 <- extract_only(prior_config[[i]])
 #      val2 <- NA
 #      print(distribution)
 #      print(paste(val1, val2))
 #    } else {
 #      stop(paste(distribution, " is not a currently supported distribution", sep = ""))
 #    }
 #
 #    # then, turn those into a data frame to match against the prior skeleton
 #    # and
 #
 # z <-  data.frame(name = c("a", "b"), value = c(1,2)) %>%
 #      spread(name, value, convert = F) %>%
 #    as.list(.)
 #  }
 #
 #  out.prior.list$p_m_center <- extract_first(prior_config[[which(names(prior_config) == "center")]])
 #  out.prior.list$p_sd_center <- extract_last(prior_config[[which(names(prior_config) == "center")]])
 #  out.prior.list$p_m_width <- extract_first(prior_config[[which(names(prior_config) == "width")]])
 #  out.prior.list$p_sd_width <- extract_last(prior_config[[which(names(prior_config) == "width")]])


    # if (stringr::str_count(prior_list[[i]], "^normal\\(") == 1) {
    #   mean <- extract_first(prior_list[[i]])
    #   sd <- extract_last(prior_list[[i]])
    #   init.name <- paste("init", names(prior_list)[i], sep = ".")
    #   assign(x = init.name, value = rlang::expr(rnorm(n = 1, mean = !!mean, sd = !!sd)), inherits = T)
    # } else if (stringr::str_count(prior_list[[i]], "^uniform\\(") == 1) {
    #   low <- extract_first(prior_list[[i]])
    #   up <- extract_last(prior_list[[i]])
    #   init.name <- paste("init", names(prior_list)[i], sep = ".")
    #   assign(x = init.name, value = rlang::expr(runif(n = 1, min = !!low, max = !!up)), inherits = T)
    # } else if (stringr::str_count(prior_list[[i]], "^exponential\\(") == 1) {
    #   rate <- extract_only(prior_list[[i]])
    #   init.name <- paste("init", names(prior_list)[i], sep = ".")
    #   assign(x = init.name, value = rlang::expr(rexp(n = 1, rate = !!rate)), inherits = T)
    # } else {
    #   message <- paste0("The distribution you've selected for parameter\n",
    #                     names(prior_list[i]), "\n",
    #                     "Isn't supported yet. Submit a feature request on github!")
    #   stop(message)
    # }

  # out.prior.list
}

