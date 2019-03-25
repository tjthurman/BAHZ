#' bahz: Bayesian Analysis of Hybrid Zones
#'
#' @description A package for fitting one-dimensional cline models on data from hybrid zones.
#'
#' @docType package
#' @name bahz-packge
#' @aliases bahz
#' @useDynLib bahz, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.4. http://mc-stan.org
#'

NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

