#' Plot a genetic cline
#'
#' Descriptions
#'
#' Details
#'
#' @importClassesFrom rstan stanfit
#'
#' @importFrom graphics "plot" "title" "points"
#'
#' @param stanfit A \code{\linkS4class{stanfit}} object holding your model results.
#'
#' @param data The dataframe with your cline data.
#'
#' @param num.out Optional, the number of x values (distances) at which to
#'   evaluate the cline. By default, does twice the length of the cline, but you
#'   can specify more or fewer. Too few may lead to a jagged-looking cline.
#'
#' @param add.obs.freqs Should the observed allele frequencies at each site be
#'   plotted? TRUE or FALSE, default is FALSE.
#'
#' @param point.col The color to use for plotting the observed allele
#'   frequencies. Default is black.
#'
#' @param ... Further graphical parameters to be passed to the base R plotting
#'   functions, in particular graphical parameters (see \code{\link[graphics]{par}}).
#'
#' @return invisible(NULL)
#'
#' @examples
#' \dontrun{
#'
#' # Default plot with the cline only
#' plot_geno_cline(yourStanfit, data)
#'
#'
#' # If you want to specify only 100 data points for plotting
#' predict_genocline(yourStanfit, data, num.out = 100)
#' }
#'
#' @export


plot_geno_cline <- function(stanfit, data, num.out = NULL, add.obs.freqs = F, point.col = "black", ...) {

  # Check arguments
  assertthat::assert_that(add.obs.freqs %in% c(T, F),
                          msg = "add.obs.freq must be True or False")

  assertthat::assert_that(point.col %in% colors(),
                          msg = paste(point.col,
                                      " is not a valid color name", sep = ""))
  # all other args get checked in predict geno cline
  cline <- bahz::predict_geno_cline(stanfit, data, num.out = num.out)

  if (add.obs.freqs == T) {
    dataframe <- data
    if (sum(c("nFocalAllele", "nTotalAlleles", "transectDist") %in% names(dataframe)) == 3) {
      assertthat::assert_that(is.numeric(dataframe$transectDist) == T,
                              msg = "transectDist column must be numeric")
      assertthat::assert_that(is.integer(dataframe$nFocalAllele) == T,
                              msg = "nFocalAllele column must be integer")
      assertthat::assert_that(is.integer(dataframe$nTotalAlleles) == T,
                              msg = "nFocalAllele column must be integer")

      # estimate allele freqs
      conv.dataframe <- dataframe %>%
        dplyr::mutate(est.p.internal = .data$nFocalAllele/.data$nTotalAlleles)


    } else if (sum(c("AA", "Aa", "aa", "transectDist") %in% names(dataframe)) == 4) {
      assertthat::assert_that(is.numeric(dataframe$transectDist) == T,
                              msg = "transectDist column must be numeric")
      assertthat::assert_that(is.integer(dataframe$AA) == T,
                              msg = "AA column must be integer")
      assertthat::assert_that(is.integer(dataframe$Aa) == T,
                              msg = "Aa column must be integer")
      assertthat::assert_that(is.integer(dataframe$aa) == T,
                              msg = "aa column must be integer")


      # estimate allele freqs
      conv.dataframe <- dataframe %>%
        dplyr::mutate(nFocalAllele = 2*.data$AA + .data$Aa,
                      nTotalAlleles = (.data$AA+.data$Aa+.data$aa)*2) %>%
        dplyr::mutate(est.p.internal = .data$nFocalAllele/.data$nTotalAlleles)

    } else {
      stop(paste("Necessary columns for calculating and plotting allele frequencies not found in data\n",
                 "Make sure it is the same data frame you used to generate the cline fit", sep = ""))
    }
  }

  graphics::plot(cline$transectDist, cline$p, type = "l", ann = F, ylim = c(0,1), ...)
  if (add.obs.freqs == T) {
    graphics::points(x = conv.dataframe$transectDist,
           y = conv.dataframe$est.p.internal,
           pch = 1,
           col = point.col)
  }
  graphics::title(...)

  return(invisible(NULL))

}
