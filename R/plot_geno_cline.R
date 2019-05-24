#' Plot a genetic cline
#'
#' Makes a simple plot to visualize modeled genetic clines. The user supplies the stanfit
#' object containing the model fit and the dataframe with the original data, and
#' may also supply a number of optional arguments to customize the plot.
#'
#' This plotting function is mostly a wrapper around
#' \code{\link{predict_geno_cline}}. For greater customization of plots, users
#' are encouraged to use \code{\link{predict_geno_cline}} to generate the x- and
#' y-coordinates for their fitted cline, and then graph those coordinates using
#' the plotting methods and packages of their choice (base plotting, lattice, or
#' ggplot2).
#'
#' @importClassesFrom rstan stanfit
#'
#' @importFrom graphics "plot" "title" "points"
#'
#' @param stanfit A \code{\linkS4class{stanfit}} object holding your model
#'   results.
#'
#' @param data The dataframe with your cline data (ideally, the same data frame
#'   that was used to generate the model fit).
#'
#' @param add.obs.freqs Should the observed allele frequencies at each site be
#'   plotted? TRUE or FALSE, default is FALSE.
#'
#' @param point.col The color to use for plotting the observed allele
#'   frequencies. Default is black.
#'
#' @param ... Further graphical parameters to be passed to the base R plotting
#'   functions to customize the plot (see \code{\link[graphics]{par}}).
#'
#' @return invisible(NULL)
#'
#' @examples
#'
#' \dontrun{
#'
#' # Default plot with the cline only
#' plot_geno_cline(yourStanfit, data)
#'
#' # Add points showing the empirical allele frequencies at
#' # each collecting site
#' plot_geno_cline(yourStanfit, data, add.obs.freq = T)
#'
#' # Some plot customization
#' # Adding axis labels, titles, and changing the
#' # colors of the points and line.
#' plot_geno_cline(yourStanfit, data, add.obs.freq = T,
#'                 main = "My cline",
#'                 xlab = "distance",
#'                 ylab = "allele frequency",
#'                 point.col = "red",
#'                 col = "blue")
#' }
#'
#' @export


plot_geno_cline <- function(stanfit, data, add.obs.freqs = F, point.col = "black", ...) {

  # Check arguments
  assertthat::assert_that(add.obs.freqs %in% c(T, F),
                          msg = "add.obs.freq must be True or False")

  assertthat::assert_that(point.col %in% colors(),
                          msg = paste(point.col,
                                      " is not a valid color name", sep = ""))
  # all other args get checked in predict geno cline
  cline <- bahz::predict_geno_cline(stanfit, data)

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
