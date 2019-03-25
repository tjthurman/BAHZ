#' Prepare your genetic cline data for loading into Stan
#'
#' Converts a dataframe containing your genetic data into the list format required for
#' input to Stan. See details below for a description of the possible input formats.
#'
#'
#' The input dataframe can be a data frame or a tibble. For the multinomial
#' model, four named columns must be present, all other columns will be
#' ignored:
#' \itemize{
#'     \item transectDist: A numeric column, giving the position along
#'     the cline/transect for each site.
#'     \item AA: The number of sampled individuals that are homozygous for the
#'     focal allele. Integer.
#'     \item Aa: The number of sampled individuals that are heterozygotes.
#'     Integer.
#'     \item aa: The number of sampled individuals that are homozygous for the
#'     non-focal allele. Integer.
#' }
#'
#' For the binomial model, the user can supply either a data frame with the four
#' columns above, or a dataframe with three named columns present (again, other
#' columns are ignored):
#'
#' \itemize{
#'     \item transectDist: A numeric column, giving the position along
#'     the cline/transect for each site.
#'     \item nFocalAllele: The allele count for the focal allele.
#'     \item nTotalAlleles: The total number of alleles sampled (twice the
#'     number of diploid individiduals).
#' }
#'
#'
#' @param dataframe A dataframe containing your cline data. See details for
#'   possible formats.
#' @param type Model type. Either "bi", for the binomial model, or "multi", for
#'   the multinomial model.
#'
#' @return Your data in Stan-ready list format.
#'
#' @seealso \code{\link{fit_geno_cline}}
#'
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' prep_geno_data(yourdata, type = "bi")
#' }
#'

prep_geno_data <- function(dataframe, type = c("bi", "multi")) {
  assertthat::assert_that(is.data.frame(dataframe) == T, msg = "dataframe must be a data frame or tibble")
  if (type == "bi") {
    if (sum(c("nFocalAllele", "nTotalAlleles", "transectDist") %in% names(dataframe)) == 3) {
      assertthat::assert_that(is.numeric(dataframe$transectDist) == T,
                              msg = "transectDist column must be numeric")
      assertthat::assert_that(is.integer(dataframe$nFocalAllele) == T,
                              msg = "nFocalAllele column must be integer")
      assertthat::assert_that(is.integer(dataframe$nTotalAlleles) == T,
                              msg = "nFocalAllele column must be integer")

      # Guess whether cline is decreasing:
      # Find row contain the first and last site on the transect
      f <- which(dataframe$transectDist == min(dataframe$transectDist))
      l <- which(dataframe$transectDist == max(dataframe$transectDist))
      # estimate allele freqs and determine if decreasing
      freq.f <- dataframe$nFocalAllele[f]/dataframe$nTotalAlleles[f]
      freq.l <- dataframe$nFocalAllele[l]/dataframe$nTotalAlleles[l]

      if (freq.f == freq.l) {
        stop("Allele frequencies equal at start and end of transect\nThey must be different for BAHZ to determine if cline is increasing or decreasing")
      }
      decrease <- as.integer(freq.f > freq.l)

      # Convert to list for Stan
      data.list <- with(dataframe, list(N = length(transectDist),
                                        nFocalAllele = nFocalAllele,
                                        nTotalAlleles = nTotalAlleles,
                                        transectDist = transectDist,
                                        decrease = decrease))

    } else if (sum(c("AA", "Aa", "aa", "transectDist") %in% names(dataframe)) == 4) {
      assertthat::assert_that(is.numeric(dataframe$transectDist) == T,
                              msg = "transectDist column must be numeric")
      assertthat::assert_that(is.integer(dataframe$AA) == T,
                              msg = "AA column must be integer")
      assertthat::assert_that(is.integer(dataframe$Aa) == T,
                              msg = "Aa column must be integer")
      assertthat::assert_that(is.integer(dataframe$aa) == T,
                              msg = "aa column must be integer")

      # Guess whether cline is decreasing:
      # Find row contain the first and last site on the transect
      f <- which(dataframe$transectDist == min(dataframe$transectDist))
      l <- which(dataframe$transectDist == max(dataframe$transectDist))
      # estimate allele freqs
      conv.dataframe <- dataframe %>%
        dplyr::mutate(nFocalAllele = 2*.data$AA + .data$Aa,
                      nTotalAlleles = (.data$AA+.data$Aa+.data$aa)*2) %>%
        dplyr::mutate(est.p.internal = .data$nFocalAllele/.data$nTotalAlleles)
      # If equal, throw a warning:
      if (conv.dataframe$est.p.internal[f] == conv.dataframe$est.p.internal[l]) {
        stop("Allele frequencies equal at start and end of transect\nThey must be different for BAHZ to determine if cline is increasing or decreasing")
      }
      decrease <- as.integer(conv.dataframe$est.p.internal[f] > conv.dataframe$est.p.internal[l])

      data.list <- with(conv.dataframe, list(N = length(transectDist),
                                        nFocalAllele = nFocalAllele,
                                        nTotalAlleles = nTotalAlleles,
                                        transectDist = transectDist,
                                        decrease = decrease))
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

      # Guess whether cline is decreasing:
      # Find row contain the first and last site on the transect
      f <- which(dataframe$transectDist == min(dataframe$transectDist))
      l <- which(dataframe$transectDist == max(dataframe$transectDist))
      # estimate allele freqs and determine if decreasing
      conv.dataframe <- dataframe %>%
        dplyr::mutate(nFocalAllele = 2*.data$AA + .data$Aa,
                      nTotalAlleles = (.data$AA+.data$Aa+.data$aa)*2) %>%
        dplyr::mutate(est.p.internal = .data$nFocalAllele/.data$nTotalAlleles)

      # If equal, throw a warning:
      if (conv.dataframe$est.p.internal[f] == conv.dataframe$est.p.internal[l]) {
        stop("Allele frequencies equal at start and end of transect\n  They must be different for BAHZ to determine if cline is increasing or decreasing")
      }

      decrease <- as.integer(conv.dataframe$est.p.internal[f] > conv.dataframe$est.p.internal[l])


      data.list <- with(dataframe, list(N = length(transectDist),
                                        genos = as.matrix(cbind(AA, Aa, aa)),
                                        transectDist = transectDist,
                                        decrease = decrease))
    } else {
      stop("Necessary data columns for multinomial model not found")
    }
  } else{
    stop("\ntype must be either\n`bi` for binomial model, or\n`multi` for multinomial model")
  }
  data.list
}
