#' Plot a cline
#'
#' Makes a simple plot to visualize modeled genetic or phenotypic clines. The user supplies the stanfit
#' object containing the model fit and the dataframe with the original data, and
#' may also supply a number of optional arguments to customize the plot.
#'
#' This plotting function is mostly a wrapper around
#' \code{\link{predict_cline}}. For greater customization of plots, users
#' are encouraged to use \code{\link{predict_cline}} to generate the x- and
#' y-coordinates for their fitted cline and confidence intervals,
#' and then graph those coordinates using
#' the plotting methods and packages of their choice (base plotting, lattice, or
#' ggplot2).
#'
#' @importClassesFrom rstan stanfit
#'
#' @importFrom graphics "plot" "title" "points"
#'
#' @importFrom grDevices "colors"
#'
#' @importFrom scales "alpha"
#'
#' @importFrom stringr "str_detect"
#'
#' @param stanfit A \code{\linkS4class{stanfit}} object holding your model
#'   results.
#'
#' @param data The dataframe with your cline data (ideally, the same data frame
#'   that was used to generate the model fit).
#'
#' @param best.fit.line  The point estimates to use for drawing the best fit
#'   line. Either "mean" for the mean of the posterior distribution of each
#'   paramer, or "median" for the median of the posterior distribution of each
#'   parameter.
#'
#' @param add.obs Should the observed allele frequencies or trait values be
#'   plotted? TRUE or FALSE, default is FALSE.
#'
#' @param cline.col the color of the cline and confidence intervals. Default is
#'   black.
#'
#' @param point.col The color to use for plotting the observed trait values.
#'   Default is black.
#'
#' @param confidence Display credible intervals around the cline? TRUE or FALSE,
#'   default FALSE.
#'
#' @param prob The probability interval to calculate around the cline. Default
#'   is .95. Numeric, between 0 and 1.
#'
#' @param method The method for calculating credible intervals. Either "ET" for
#'   equal-tail probability intervals, or "HPDI" for highest posterior density
#'   intervals. Default is "HPDI".
#'
#' @param clear.cache Clear the cache of saved results to ensure recalculation predicted cline?
#' TRUE or FALSE, default FALSE.
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
#' plot_cline(yourStanfit, data)
#'
#' # Add points showing the observed raw data
#' plot_cline(yourStanfit, data, add.obs = T)
#'
#' # Add credible intervals around the cline
#' plot_cline(yourStanfit, data, confidence = T)
#'
#' # Some plot customization
#' # Adding axis labels, titles, and changing the
#' # colors of the points and line.
#' plot_cline(yourStanfit, data, add.obs = T,
#'                 main = "My cline",
#'                 xlab = "distance",
#'                 ylab = "trait value",
#'                 point.col = "red",
#'                 col = "blue")
#' }
#'
#' @export
#'

plot_cline <- function(stanfit, data, best.fit.line = "mean",
                       add.obs = F, confidence = F,
                       prob = 0.95,
                       method = "HPDI",
                       cline.col = "black", point.col = "black",
                       clear.cache = F, ...) {
  # Check arguments
  assertthat::assert_that(class(stanfit)[1] == "stanfit",
                          msg = "Model object from which to plot must be of class stanfit")
  assertthat::assert_that(is.data.frame(data),
                          msg = paste("Input data must be a data frame",
                                      "Make sure it is the same data frame you used to generate the cline fit",
                                      sep = "\n"))
  match.arg(best.fit.line, choices = c("mean", "median"), several.ok = F)
  assertthat::assert_that("transectDist" %in% colnames(data),
                          msg = paste("Input data frame does not contain a transectDist column",
                                      "Make sure it is the same data frame you used to generate the cline fit",
                                      sep = "\n"))
  assertthat::assert_that(is.numeric(data$transectDist),
                          msg = paste("transectDist column in input data must be numeric",
                                      "Make sure it is the same data frame you used to generate the cline fit",
                                      sep = "\n"))
  assertthat::assert_that(is.logical(add.obs),
                          msg = "add.obs must be True or False")
  assertthat::assert_that(cline.col %in% colors(),
                          msg = paste(point.col,
                                      " is not a valid color name", sep = ""))
  assertthat::assert_that(point.col %in% colors(),
                          msg = paste(point.col,
                                      " is not a valid color name", sep = ""))
  assertthat::assert_that((method %in% c("HPDI", "ET")) == T,
                          msg = "method must be either 'HPDI' or 'ET")

  # Check supplied extra graphical parameters
  extra.args <- list(...)
  reserved.par <- c("ann", "border", "col", "type", "ylim")
  if (sum(reserved.par %in% names(extra.args)) > 0) {
    errs <- reserved.par[which(reserved.par %in% names(extra.args))]
    if ("col" %in% errs) {
      msg <- paste("Some graphical parameters supplied in ...",
                   "have default values in bahz and cannot be overridden.",
                   "\n",
                   "Remove these arguments:",
                   toString(errs),
                   "\n",
                   "Use cline.col and point.col to specify colors,",
                   "don't include col as an additional argument.", sep = "\n")
    } else {
      msg <- paste("Some graphical parameters supplied in ...",
                   "have default values in bahz and cannot be overridden:",
                   "\n",
                   "Remove these arguments:",
                   toString(errs),
                   sep = "\n")
    }
    stop(msg)
  }

  # Figure out if a pheno cline or a geno cline
  if (stringr::str_detect(stanfit@model_name, "pheno")) {
    phenotypic <- T
  } else {
    phenotypic <- F
  }



  # Check to see if the input data frame matches in the stanfit object
  # in terms of rows of raw data
  # as in the input data frame. Give a warning if there's not.
  data.rows <- dim(data)[1]
  if (phenotypic) {
    sf.rows <- stanfit@par_dims$y_rep
  } else {
    sf.rows <- stanfit@par_dims$p
  }
  if (data.rows != sf.rows) {
    warning(paste("\n",
                  "Your stanfit object and dataframe contain data with different numbers of rows\n",
                  "The stanfit object may have been generated from a different data frame\n",
                  "stanfit rows = ",
                  sf.rows,
                  "\n",
                  "data rows = ",
                  data.rows, sep = ""))
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
  cline <- bahz::predict_cline(stanfit, distance = xrange,
                               confidence = confidence,
                               prob = prob, method = method,
                               clear.cache = clear.cache)

  # If adding the observed data, check that the decessary rows are there.


  if (add.obs) {
    if (phenotypic) {
      # if a phenotypic cline.
      if (sum(c("transectDist", "traitValue") %in% names(data)) == 2) {
        assertthat::assert_that(is.numeric(data$transectDist) == T,
                                msg = "transectDist column must be numeric")
        assertthat::assert_that(is.numeric(data$traitValue) == T,
                                msg = "traitValue column must be numeric")

        conv.dataframe <- data %>%
          dplyr::arrange(.data$transectDist)
      }
      else {# If the data columns aren't present
        stop(
          paste(
            "Necessary columns for plotting observed trait values not found in data\n",
            "Make sure it is the same data frame you used to generate the cline fit",
            sep = ""
          )
        )
      }
    } else {
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
          dplyr::mutate(est.p.internal = .data$nFocalAllele / .data$nTotalAlleles)
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
          dplyr::mutate(
            nFocalAllele = 2 * .data$AA + .data$Aa,
            nTotalAlleles = (.data$AA + .data$Aa + .data$aa) *
              2
          ) %>%
          dplyr::mutate(est.p.internal = .data$nFocalAllele / .data$nTotalAlleles)

      } else {
        stop(
          paste(
            "Necessary columns for calculating and plotting allele frequencies not found in data\n",
            "Make sure it is the same data frame you used to generate the cline fit",
            sep = ""
          )
        )
      }
    }
  }

  # Now, plot it all
  if (phenotypic) {
    ylims <- c(min(data$traitValue), max(data$traitValue))
  } else {
    ylims <- c(0,1)
  }
  if (best.fit.line == "median") {
    graphics::plot(cline$transectDist, cline$p_median, type = "l", ann = F, col = cline.col, ylim = ylims, ...)
  } else {
    graphics::plot(cline$transectDist, cline$p_mean, type = "l", ann = F, col = cline.col, ylim = ylims, ...)

  }
  if (confidence) {
    graphics::polygon(x = c(cline$transectDist, rev(cline$transectDist)),
                      y = c(cline[,4], rev(cline[,5])), border = NA,
                      col = scales::alpha(cline.col, 0.2), ...)
    if (best.fit.line == "median") {
      graphics::lines(cline$transectDist, cline$p_median, col = cline.col, ...)
    } else {
      graphics::lines(cline$transectDist, cline$p_mean, col = cline.col, ...)
    }
  }

  if (add.obs == T) {
    if (phenotypic) {
      graphics::points(x = conv.dataframe$transectDist,
                       y = conv.dataframe$traitValue,
                       col = point.col, ...)
    }
    else {
      graphics::points(x = conv.dataframe$transectDist,
                       y = conv.dataframe$est.p.internal,
                       col = point.col, ...)
    }
  }
  graphics::title(...)

  return(invisible(NULL))
}

