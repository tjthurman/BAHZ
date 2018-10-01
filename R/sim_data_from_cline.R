#' Simulate genetic data from a cline
#'
#' This function generates a dataframe with simulated genotypic data sampled
#' from a genetic cline. The sampling sites, number of individuals, level of
#' inbreeding, and cline parameters are supplied by the user. Cline parameters
#' are flexible, and can model both sigmoid clines and stepped clines with
#' introgresison tails.
#'
#' DETAILS TO BE ADDED. SHOULD PUT EQUATIONS IN \code{\link{general_cline_eqn}} AND LINK TO IT.
#'
#' @importFrom stats "rmultinom"
#'
#' @param transect_distances The distances along the transect for the simulated
#'   sampling sites. A numeric vector.
#' @param n_ind The number of diploid individuals sampled at each site. Either a
#'   single numeric value (for constant sampling), or a numeric vector equal in
#'   length to \code{transect_distances}.
#' @param Fis The inbreeding coefficient, Fis, for each site. Must be between 0
#'   and 1 (inclusive). Either a single numeric value (for constant inbreeding),
#'   or a numeric vector equal in length to \code{transect_distances}.
#' @param center The location of the cline center, in the same distance units as
#'   \code{transect_distances}. Numeric, must be greater than 0.
#' @param width The width of the cline, in the same distance units as
#'   \code{transect_distances}. Numeric, must be greater than 0.
#' @param pmin,pmax Optional. The minimum and maximum allele frequency values
#'   in the tails of the cline. Default values are \code{0} and \code{1}, respectively.
#'   Must be between 0 and 1 (inclusive). Numeric.
#' @param deltaL,tauL Optional delta and tau parameters which describe the left
#'   exponential tail. Must supply both to generate a tail. Default is \code{NULL} (no
#'   tails). Numeric. tauL must be between 0 and 1 (inclusive).
#' @param deltaR,tauR Optional delta and tau parameters which describe the right
#'   exponential tail. Must supply both to generate a tail. Default is \code{NULL} (no
#'   tails). Numeric. tauR must be between 0 and 1 (inclusive).
#'
#' @importFrom magrittr "%>%"
#'
#' @return A data frame of simulated genetic data sampled from the cline. Columns are:
#'     \itemize{
#'     \item site: The site numbers, given sequentially starting at 1.
#'     \item transect.Dist: The distance along the cline for each site.
#'     \item cline.p: The expected allele frequency for each site, given its position on the cline.
#'     \item cline.f: The expected coefficient of inbreeding for each site.
#'     \item AA, AA, aa: The simulated number of homozygotes and heterozygotes for each site.
#'     \item N: The number of individuals sampled for each site.
#'     \item emp.p: The observed allele frequency for each site (includes sampling error).
#'     \item emp.f: The observed Fis for each site (includes sampling error).
#' }
#'
#' @export
#'
#' @examples
#' # Simulate genotype data from a cline with center at 100, width of 30.
#' # Sites are 20 units apart, from 0 to 200.
#' # 20 individuals are sampled at each site.
#' # Inbreeding is constant at Fis = 0.1.
#'
#' set.seed(123)
#' sim_data_from_cline(transect_distance =seq(0,200,20), n_ind = 20,
#'                     Fis = 0.1, decrease = T,
#'                     center = 100, width = 30)
#'

sim_data_from_cline <- function(transect_distances, n_ind,
                                Fis, decrease,
                                center, width,
                                pmin = 0, pmax = 1,
                                deltaL = NULL, tauL = NULL,
                                deltaR = NULL, tauR = NULL) {


  # Check the sampling and inbreeding options
  for (vec.arg in alist(n_ind, Fis)) {
    assertthat::assert_that(is.vector(eval(vec.arg)) == T,
                msg = paste(vec.arg, "must be a vector", sep = " "))
    assertthat::assert_that(is.numeric(eval(vec.arg)) == T,
                msg = paste(vec.arg, "must be numeric", sep = " "))
    assertthat::assert_that((length(eval(vec.arg)) %in% c(1, length(transect_distances))) == T,
               msg = paste(vec.arg, " must be either be of length 1, for constant ", vec.arg,
                           ", or must match the length of transect_distances (",
                           length(transect_distances), sep = ""))
  }
  assertthat::assert_that(min(Fis) >=0, msg = "Fis values cannot be less than 0")
  assertthat::assert_that(min(Fis) <=1, msg = "Fis values cannot be greater than 1")
  # All other args will get checked in the cline equation.

  # Get number of sites from the vector of transect data.
  sites <- length(transect_distances)
  # Get the vector of f values for each site
  if (length(n_ind) == 1) {
    Ns <- rep(n_ind, times = sites)
  } else {
      Ns <- n_ind
  }
  if (length(Fis) == 1) {
    fs <- rep(Fis, times = sites)
  } else {
    fs <- Fis
  }

  # Make the empty results data frame
  fk.dt <- data.frame(site = 1:sites,
                      transect.Dist = transect_distances,
                      cline.p = rep(NA, times = sites),
                      cline.f = fs,
                      AA = rep(NA, times = sites),
                      Aa = rep(NA, times = sites),
                      aa = rep(NA, times = sites),
                      N = Ns)

  # Then add the simulated genotypes to each row
  for (row in 1:sites) {
    fk.dt$cline.p[row] <- general_cline_eqn(transectDist = fk.dt$transect.Dist[row], center = center, width = width,
                                            pmin = pmin, pmax = pmax, deltaL = deltaL, deltaR = deltaR, tauL = tauL,
                                            tauR = tauR, decrease = decrease)
    AA <- fk.dt$cline.p[row]^2 + fk.dt$cline.f[row]*fk.dt$cline.p[row]*(1-fk.dt$cline.p[row])
    Aa<- 2*fk.dt$cline.p[row]*(1-fk.dt$cline.p[row])*(1-fk.dt$cline.f[row])
    aa <- (1-fk.dt$cline.p[row])^2 +fk.dt$cline.f[row]*fk.dt$cline.p[row]*(1-fk.dt$cline.p[row])
    genotypes <- rowSums(stats::rmultinom(n = 1, size = fk.dt$N[row], prob = c(AA, Aa, aa)))
    fk.dt$AA[row] <- genotypes[1]
    fk.dt$Aa[row] <- genotypes[2]
    fk.dt$aa[row] <- genotypes[3]
    fk.dt$N[row] <- sum(genotypes)
  }


# Calculate empirical p and Fis value from the simulated data
fk.dt <- fk.dt %>%
  dplyr::mutate(emp.p = (2*AA + Aa)/(2*N)) %>%
  dplyr::mutate(Hexp = 2*emp.p*(1-emp.p),
         Hobs = Aa/N) %>%
  dplyr::mutate(Fis = (Hexp - Hobs)/Hexp) %>%
  correct_fis(.) %>%
  dplyr::rename(emp.f = Fis) %>%
  dplyr::select(-Hexp, -Hobs)

  # Do some rounding
  fk.dt$cline.p <- round(fk.dt$cline.p, digits = 3)
  fk.dt$emp.p <- round(fk.dt$emp.p, digits = 3)
  fk.dt$emp.f <- round(fk.dt$emp.f, digits = 3)
  fk.dt
}
