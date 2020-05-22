
# Load packages -----------------------------------------------------------
rm(list = ls())
library(devtools)
detach("package:bahz", unload=TRUE)
install_github("tjthurman/BAHZ")
devtools::install_github("tjthurman/BAHZ@feature/dongyiyi_install_fixes")
install.packages("~/Documents/Work/PhD:McGill/Projects/BAHZ/", repos = NULL, type = "source")
library(bahz)
library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(loo)
# Getting the minimal pre-compiled model to run -------------------
# Generate a dataset
set.seed(22)
data <- sim_geno_cline(transect_distances = seq(-300,300,20), n_ind = 40, Fis = .9,
                       decrease = F, center = 10, width = 35, pmin = 0.08, pmax = .95, deltaR = 12, tauR = 0.25)

plot(x = data$transectDist, y = data$emp.p)
lines(x = data$transectDist, y = data$cline.p)


library(bahz)

prep_prior_list("~/Desktop/geno_priors.yaml")


# Fit the model

none_bi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                           type = "bi", tails = "none")
none_multi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                              type = "multi", tails = "none")

left_bi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                               type = "bi", tails = "left")
left_multi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                           type = "multi", tails = "left")

right_bi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                           type = "bi", tails = "right")
right_multi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                            type = "multi", tails = "right")

mirror_bi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                            type = "bi", tails = "mirror")
mirror_multi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                             type = "multi", tails = "mirror")

ind_bi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                             type = "bi", tails = "ind")
ind_multi <- fit_geno_cline2(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                          type = "multi", tails = "ind")



rethinking::compare(none_bi, left_bi, right_bi, mirror_bi, ind_bi)
cline_summary(none_bi)
cline_summary(none_multi)
cline_summary(left_bi)
cline_summary(left_multi)
cline_summary(right_bi)
cline_summary(right_multi)
cline_summary(mirror_bi)
cline_summary(mirror_multi)
cline_summary(ind_bi)
cline_summary(ind_multi)


geno_fit_bi@stanmodel

geno_fit_bi@sim$samples
library(bayesplot)
color_scheme_set("viridis")
mcmc_trace(as.array(geno_fit_bi), pars = c("deltaL", "deltaR"))
mcmc_trace(as.array(geno_fit_bi), pars = c("deltaL"))
mcmc_trace(as.array(geno_fit_bi2), pars = c("center"))


cline_summary(geno_fit_bi)
geno_fit_bi2 <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                              type = "bi", tails = "left")

pairs(geno_fit_bi)

geno_fit_bi@inits


log(2)
cline_summary(geno_fit_bi)
cline_summary(geno_fit_bi2)


geno_fit_multi <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                           type = "multi", tails = "none")

cline_summary(geno_fit_bi)

stanfit <- geno_fit_bi




# Will stan run with no priors? -------------------------------------------

# Yes, stan will run with no priors! I think the default, then, is that it
# puts uniform -inf, Inf priors, which against the whole point of using stan,
# really. Still, good to test.


# Re-scaling of width? ----------------------------------------------------

# Is it possible to rescale the width parameter as well, using 4/width

# Testing on how to convert from normal distribution on width scale
# to normal distribution on w scale.

# After a lot of testing, the answer seems to be no.
# It would be nice to re-scale width to ease the starting and allow it
# to start from default values. But, I'm nt sure there's any way to convert
# the user-supplied prior to a sensible prior on the 4/width scale.


