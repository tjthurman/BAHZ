#' Plot a genetic cline
#'
#' Makes a simple plot to visualize modeled genetic clines. The user supplies the stanfit
#' object containing the model fit and the dataframe with the original data, and
#' may also supply a number of optional arguments to customize the plot.
#'
#' This plotting function is mostly a wrapper around
#' \code{\link{predict_cline}}. For greater customization of plots, users
#' are encouraged to use \code{\link{predict_cline}} to generate the x- and
#' y-coordinates for their fitted cline, and then graph those coordinates using
#' the plotting methods and packages of their choice (base plotting, lattice, or
#' ggplot2).
#'
#' @importClassesFrom rstan stanfit
#'
#' @importFrom graphics "plot" "title" "points"
#'
#' @importFrom grDevices "colors"
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
  assertthat::assert_that(class(stanfit)[1] == "stanfit",
                          msg = "Model object from which to plot must be of class stanfit")
  assertthat::assert_that(is.data.frame(data),
                          msg = paste("Input data must be a data frame",
                                      "Make sure it is the same data frame you used to generate the cline fit",
                                      sep = "\n"))
  assertthat::assert_that("transectDist" %in% colnames(data),
                          msg = paste("Input data frame does not contain a transectDist column",
                                      "Make sure it is the same data frame you used to generate the cline fit",
                                      sep = "\n"))
  assertthat::assert_that(is.numeric(data$transectDist),
                          msg = paste("transectDist column in input data must be numeric",
                                      "Make sure it is the same data frame you used to generate the cline fit",
                                      sep = "\n"))
  assertthat::assert_that(add.obs.freqs %in% c(T, F),
                          msg = "add.obs.freq must be True or False")
  assertthat::assert_that(point.col %in% colors(),
                          msg = paste(point.col,
                                      " is not a valid color name", sep = ""))

  # Check to see if there are the same number of sites in the stanfit object
  # as in the inout data frame. Give a warning if there's not.
  data.sites <- dim(data)[1]
  sf.sites <- stanfit@par_dims$p
  if (data.sites != sf.sites) {
    warning(paste("\n",
                  "Your stanfit object and dataframe contain data from different numbers of sites\n",
                  "The stanfit object may have been generated from a different data frame\n",
                  "stanfit sites = ",
                  sf.sites,
                  "\n",
                  "data sites = ",
                  data.sites, sep = ""))
  }

  # Figure out starting and ending x values
  start.x <- as.integer(min(data$transectDist))
  end.x <- as.integer(max(data$transectDist))
  range <- end.x - start.x
  start.pad <- as.integer(start.x - (0.02*(range)))
  end.pad <- as.integer(end.x + (0.02*(range)))

  xrange <- seq(from = start.pad, to = end.pad,
                length.out = range*2)

  # Generate the cline values to plot
  cline <- bahz::predict_cline(stanfit, distance = xrange)

  # If addinf the observed allele frequencies, calculate those.
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

  # Now, plot it all
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
