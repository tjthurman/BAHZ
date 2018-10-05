#' Prepare your cline data for loading into Stan
#'
#' DESCRIPTION TO BE WRITTEN
#'
#' DETAILS TO BE ADDED. SPECIFICALLY, REQUIRED FILE FORMAT.
#'
#' @param dataframe A dataframe containing your cline data. See details for
#'   possible formats.
#' @param type Model type. Either "bi", for the binomial model, or "multi", for
#'   the multinomial model.
#'
#' @return Your data, in the list format appropriate for the cline model you
#'   chose.
#'
#' @export
#'
#' @examples
#' # TO BE ADDED
#'



load_cline_data <- function(dataframe, type) {
  if (type == "bi") {
    if (sum(c("nFocalAllele", "nTotalAlleles", "transectDist") %in% names(dataframe)) == 3) {
      assertthat::assert_that(is.numeric(dataframe$transectDist) == T,
                              msg = "transectDist column must be numeric")
      assertthat::assert_that(is.integer(dataframe$nFocalAllele) == T,
                              msg = "nFocalAllele column must be integer")
      assertthat::assert_that(is.integer(dataframe$nTotalAlleles) == T,
                              msg = "nFocalAllele column must be integer")

      data.list <- with(dataframe, list(N = length(transectDist),
                                        nFocalAllele = nFocalAllele,
                                        nTotalAlleles = nTotalAlleles,
                                        transectDist = transectDist))
    } else if (sum(c("AA", "Aa", "aa", "transectDist") %in% names(dataframe)) == 4) {
      assertthat::assert_that(is.numeric(dataframe$transectDist) == T,
                              msg = "transectDist column must be numeric")
      assertthat::assert_that(is.integer(dataframe$AA) == T,
                              msg = "AA column must be integer")
      assertthat::assert_that(is.integer(dataframe$Aa) == T,
                              msg = "Aa column must be integer")
      assertthat::assert_that(is.integer(dataframe$aa) == T,
                              msg = "aa column must be integer")
      data.list <- with(dataframe, list(N = length(transectDist),
                                        nFocalAllele = 2*AA + Aa,
                                        nTotalAlleles = 2*N,
                                        transectDist = transectDist))
    } else {
      stop("Necessary data columns for binomial model not found")
    }
  } else if (type == "multi") {
    if (sum(c("AA", "Aa", "aa", "transectDist") %in% names(dataframe)) == 4) {
      assertthat::assert_that(is.numeric(dataframe$transectDist) == T,
                              msg = "transectDist column must be numeric")
      assertthat::assert_that(is.integer(dataframe$AA) == T,
                              msg = "AA column must be integer")
      assertthat::assert_that(is.integer(dataframe$Aa) == T,
                              msg = "Aa column must be integer")
      assertthat::assert_that(is.integer(dataframe$aa) == T,
                              msg = "aa column must be integer")
      data.list <- with(dataframe, list(N = length(transectDist),
                                        genos = as.matrix(cbind(AA, Aa, aa)),
                                        transectDist = transectDist))
    } else {
      stop("Necessary data columns for multinomial model not found")
    }
  } else{
    stop("\ntype must be either\n`bi` for binomial model, or\n`multi` for multinomial model")
  }
  data.list
}
