#' Prepare your phenotypic cline data for loading into Stan
#'
#' Converts a dataframe containing your phenotypic data into the list format required for
#' input to Stan. See details below for a description of the possible input formats.
#'
#'
#' The input dataframe can be a data frame or a tibble. Each row should contain
#' the information for a single sampled individual. Two named
#' columns must be present, all other columns will be ignored:
#' \itemize{
#'     \item transectDist: A numeric column, giving the position along the
#'     cline/transect for each individual.
#'     \item traitValue: A numeric column
#'     giving the phenotypic trait value for each individual.
#' }
#'
#'
#' @param dataframe A dataframe containing your phenotypic cline data. See
#'   details for the necessary format.
#'
#' @return Your data in Stan-ready list format.
#'
#' @seealso \code{\link{fit_pheno_cline}}
#'
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' prep_pheno_data(yourdata)
#' }
#'

prep_pheno_data <- function(dataframe) {
  # Type checking
  assertthat::assert_that(is.data.frame(dataframe) == T, msg = "dataframe must be a data frame or tibble")

  # check that the necessary columns are present and make the input date
  if (sum(c("transectDist", "traitValue") %in% names(dataframe)) == 2) {
      assertthat::assert_that(is.numeric(dataframe$transectDist) == T,
                              msg = "transectDist column must be numeric")
      assertthat::assert_that(is.numeric(dataframe$traitValue) == T,
                            msg = "traitValue column must be numeric")

      dataframe <- dataframe %>%
        dplyr::arrange(.data$transectDist)

      # Guess whether cline is decreasing:
      # Find row containing the first and last site on the transect
      site_means <- dataframe %>%
        dplyr::group_by(.data$transectDist) %>%
        dplyr::summarize(site_means = mean(.data$traitValue),
                         n_per_site = length(.data$traitValue))

      f <- which(site_means$transectDist == min(dataframe$transectDist))
      l <- which(site_means$transectDist == max(dataframe$transectDist))

      freq.f <- site_means$site_means[f]
      freq.l <-  site_means$site_means[l]

      if (freq.f == freq.l) {
        stop("Average phenotypes are equal at the start and end of transect\nThey must be different for BAHZ to determine if cline is increasing or decreasing")
      }
      decrease <- as.integer(freq.f > freq.l)

      # Convert to list for Stan
      data.list <- with(dataframe, list(N = length(traitValue),
                                        K = length(unique(transectDist)),
                                        pheno = traitValue,
                                        s = site_means$n_per_site,
                                        transectDist = unique(transectDist),
                                        decrease = decrease))
    }

  # If the data columns aren't present
  else {
      missing <- c("transectDist", "traitValue")[which((c("transectDist", "traitValue") %in% names(dataframe)) == F)]

      stop(paste("\n",
                 "Necessary data columns for phenotypic cline model not found",
                 " ",
                 "Missing columns are:",
                 paste(missing, collapse = "\n"),
                 " ",
                 "See ?prep_pheno_data for more info",
                 "on the proper formatting of",
                 "input data", sep = "\n"))
    }
  data.list
}
