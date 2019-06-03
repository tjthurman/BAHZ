#' Simulate phenotypic data from a cline
#'
#' This function generates a dataframe with simulated trait data sampled from a
#' phenotypic cline. The sampling sites, number of individuals, phenotypic
#' variance (sigma), and cline parameters are supplied by the user.
#'
#' This function calls \code{\link{general_cline_eqn}} for each sampled point along the
#' cline to generate the expected mean phenotype at that point.
#'
#' For each site the specified number of individuals are simulated. Trait
#' values are drawn from a normal distribution, with mean equal to the expected mean
#' phenotype and standard devitation equal to the provided sigma.
#'
#' @importFrom magrittr "%>%"
#'
#' @param transect_distances The distances along the transect for the simulated
#'   sampling sites. A numeric vector.
#' @param n_ind The number of individuals sampled at each site. Either a
#'   single numeric value (for constant sampling), or a numeric vector equal in
#'   length to \code{transect_distances}.
#' @param sigma The simulated phenotypic variance, given in terms of standard
#'   deviation for a normal distribution. Must be numeric and > 0.  Either a
#'   single numeric value (for constant variance), or a numeric vector equal
#'   in length to \code{transect_distances}.
#' @param decrease Is the cline decreasing in frequency? \code{TRUE} or
#'   \code{FALSE}.
#' @param center The location of the cline center, in the same distance units as
#'   \code{transect_distances}. Numeric, must be greater than 0.
#' @param width The width of the cline, in the same distance units as
#'   \code{transect_distances}. Numeric, must be greater than 0.
#' @param pmin,pmax The minimum and maximum mean phenotypic trait value
#'   in the tails of the cline. Numeric.
#'
#' @return A data frame of simulated phenotypic data sampled from the cline. Each
#'   row representes a single individual. Columns are:
#'     \itemize{
#'     \item site: The site numbers, given sequentially starting at 1.
#'     \item transectDist: The distance along the cline for each individual
#'     \item traitValue: The simulated trait value for each individual.
#' }
#'
#' @export
#'
#' @examples
#' # Simulate phenotype data from a decreasing cline
#' # with center at 100, width of 30, pmin of 10 and pmax of 40.
#' # Sites are 20 units apart, from 0 to 200.
#' # 20 individuals are sampled at each site.
#' # Variance is constant at sigma = 6.
#'
#' set.seed(123)
#' sim_pheno_cline(transect_distance = seq(0,200,20),
#'                n_ind = 20, sigma = 6,
#'                decrease = TRUE,
#'                center = 100, width = 30,
#'                pmin = 10, pmax = 40)
#'
#' # Simulate phenotype data from an increasing cline
#' # with center at 272, width of 91.
#' # The minimum and maximum mean phenotype values
#' # are 32 and 49, respectively.
#'
#' # Sites are 13 units apart, from 162 to 312.
#' # At each site, the number of individuals sampled
#' # is drawn from a random normal distribution with
#' # mean = 25 and sd = 5
#' # variance is constant at sigma = 4.67.
#'
#' set.seed(123)
#' ind_sampling <- as.integer(rnorm(length(seq(162,312,12)), 25, 5))
#' sim_pheno_cline(transect_distance = seq(162,312,12),
#'                n_ind = ind_sampling,
#'                sigma = 4.67, decrease = FALSE,
#'                center = 272, width = 91,
#'                pmin = 32, pmax = 49)
#'

sim_pheno_cline <- function(transect_distances, n_ind,
                           sigma, decrease,
                           center, width,
                           pmin, pmax) {


  # Check the sampling and inbreeding options
  for (vec.arg in alist(n_ind, sigma)) {
    assertthat::assert_that(is.vector(eval(vec.arg)) == T,
                            msg = paste(vec.arg, "must be a vector", sep = " "))
    assertthat::assert_that(is.numeric(eval(vec.arg)) == T,
                            msg = paste(vec.arg, "must be numeric", sep = " "))
    assertthat::assert_that((length(eval(vec.arg)) %in% c(1, length(transect_distances))) == T,
                            msg = paste(vec.arg, " must be either be of length 1, for constant ", vec.arg,
                                        ", or must match the length of transect_distances (",
                                        length(transect_distances), sep = ""))
  }
  assertthat::assert_that(min(sigma) >=0, msg = "sigma cannot be less than 0")
  assertthat::assert_that(min(n_ind) >=1, msg = "n_ind values cannot be less than 1")
  assertthat::assert_that(pmin < pmax, msg = "pmin must be less than pmax")
  # All other args will get checked in the cline equation.

  # Get number of sites from the vector of transect data.
  sites <- length(unique(transect_distances))
  # Get the vector of individuals and sigma values for each site
  if (length(n_ind) == 1) {
    Ns <- as.integer(rep(n_ind, times = sites))
  } else {
    Ns <- as.integer(n_ind)
  }
  if (length(sigma) == 1) {
    sig <- rep(sigma, times = sites)
  } else {
    sig <- sigma
  }

  # Get the data for each site
  sim_site <- data.frame(site = 1:sites,
                         transectDist = transect_distances,
                         inds = Ns,
                         sigma = sig)
  sim_site$pred.mean <- general_cline_eqn(transectDist = sim_site$transectDist,
                                          center = center,
                                          width = width,
                                          pmin = pmin,
                                          pmax = pmax,
                                          decrease = decrease)

  # Then loop over number of sites to generate individuals
  ## NEED TO FINISH THIS

  sim_inds <- NULL

  for (site in 1:dim(sim_site)[1])
    {
    one.site <- data.frame(site = sim_site$site[site],
                           transectDist = sim_site$transectDist[site],
                           traitValue = rnorm(n = sim_site$inds[site],
                                              mean = sim_site$pred.mean[site],
                                              sd = sim_site$sigma[site]))
    sim_inds <- rbind(sim_inds, one.site)
    }

  sim_inds
}
