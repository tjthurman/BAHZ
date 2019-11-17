#' Plot a phenotypic cline
#'
#' Makes a simple plot to visualize modeled phenotypic clines. The user supplies the stanfit
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
#' @param add.obs.pheno Should the observed trait values for each individual be
#'   plotted? TRUE or FALSE, default is FALSE.
#'
#' @param point.col The color to use for plotting the observed trait
#'   values. Default is black.
#'
#' @param confidence Display credible intervals around the cline? TRUE or FALSE, default FALSE.
#'
#' @param prob The probability interval to calculate for the cline. Default is .95. Numeric,
#'   between 0 and 1.
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
#' plot_pheno_cline(yourStanfit, data)
#'
#' # Add points showing the observed trait values for each individual
#' plot_pheno_cline(yourStanfit, data, add.obs.pheno = T)
#'
#' # Some plot customization
#' # Adding axis labels, titles, and changing the
#' # colors of the points and line.
#' plot_pheno_cline(yourStanfit, data, add.obs.pheno = T,
#'                 main = "My cline",
#'                 xlab = "distance",
#'                 ylab = "trait value",
#'                 point.col = "red",
#'                 col = "blue")
#' }
#'
#' @export


plot_pheno_cline <- function(stanfit, data, add.obs.pheno = F, point.col = "black", ...) {

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
  assertthat::assert_that(add.obs.pheno %in% c(T, F),
                          msg = "add.obs.pheno must be True or False")
  assertthat::assert_that(point.col %in% colors(),
                          msg = paste(point.col,
                                      " is not a valid color name", sep = ""))
  assertthat::assert_that(is.numeric(prob) == T, msg = "prob must be numeric")
  assertthat::assert_that(length(prob) == 1, msg = "prob must be of length 1")
  assertthat::assert_that(prob <= 1, msg = "prob must be between 0 and 1")
  assertthat::assert_that(prob > 0, msg = "prob must be between 0 and 1")
  assertthat::assert_that(is.logical(confidence) == T, msg = "confidence must be either TRUE or FALSE")


  # Check to see if there are the same number of individuals in the stanfit object
  # as in the input data frame. Give a warning if there's not.
  data.inds <- dim(data)[1]
  sf.inds <- stanfit@par_dims$y_rep
  if (data.inds != sf.inds) {
    warning(paste("\n",
                  "Your stanfit object and dataframe contain data from different numbers of individuals\n",
                  "The stanfit object may have been generated from a different data frame\n",
                  "stanfit individuals = ",
                  sf.inds,
                  "\n",
                  "data individuals = ",
                  data.inds, sep = ""))
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
  if (add.obs.pheno == T) {
    if (sum(c("transectDist", "traitValue") %in% names(data)) == 2) {
      assertthat::assert_that(is.numeric(data$transectDist) == T,
                              msg = "transectDist column must be numeric")
      assertthat::assert_that(is.numeric(data$traitValue) == T,
                              msg = "traitValue column must be numeric")

      dataframe <- data %>%
        dplyr::arrange(.data$transectDist)

    }

  # If the data columns aren't present
  else {
    stop(paste("Necessary columns for plotting observed trait values not found in data\n",
               "Make sure it is the same data frame you used to generate the cline fit", sep = ""))
    }
  }

  # Now, plot it all
  graphics::plot(cline$transectDist, cline$p, type = "l", ann = F, ylim = c(min(data$traitValue), max(data$traitValue)), ...)
  if (add.obs.pheno == T) {
    graphics::points(x = dataframe$transectDist,
                     y = dataframe$traitValue,
                     pch = 1,
                     col = point.col)
  }
  graphics::title(...)

  return(invisible(NULL))

}
