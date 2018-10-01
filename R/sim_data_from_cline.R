#' Simulate genetic data from a cline
#'
#' DESCRIPTION TO BE ADDED
#'
#' DETAILS TO BE ADDED. LINK TO THE GENERAL CLINE EQUATION FUNCTION, MAYBE MOVE
#' THE CLINE PARAMETER DEFINITIONS THERE.
#'
#' @importFrom stats "rmultinom"
#'
#' @param transect_distances The distances along the transect for the sampling
#'   sites to be simulated. A numeric vector.
#' @param n_ind The number of diploid individuals sampled at each site. Either a
#'   single numeric value (for constant sampling), or a numeric vector equal in
#'   length to \code{transect_distances}.
#' @param Fis The inbreeding coefficient, Fis, for each site. Must be between 0
#'   and 1 (inclusive). Either a single numeric value (for constant inbreeding),
#'   or a numeric vector equal in length to \code{transect_distances}.
#'
#' @inheritParams general_cline_eqn
#'
#' @importFrom magrittr "%>%"
#'
#' @return TO ADD
#'
#' @export
#'
#' @examples
#' # to be added
#'
#'

sim_data_from_cline <- function(transect_distances, n_ind,
                                Fis, decrease = F,
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
                      transectDist = transect_distances,
                      cline.p = rep(NA, times = sites),
                      cline.f = fs,
                      AA = rep(NA, times = sites),
                      Aa = rep(NA, times = sites),
                      aa = rep(NA, times = sites),
                      N = Ns)

  # Then add the simulated genotypes to each row
  for (row in 1:sites) {
    fk.dt$cline.p[row] <- general_cline_eqn(transectDist = fk.dt$transectDist[row], center = center, width = width,
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
