
# Load packages -----------------------------------------------------------
rm(list = ls())
library(devtools)
detach("package:bahz", unload=TRUE)
install_github("tjthurman/BAHZ", ref = "master", auth_token = "c037ee28d7031c8d2cfb91dbefc686409bf05e68")
install.packages("~/Documents/Work/PhD:McGill/Projects/BAHZ/", repos = NULL, type = "source")
library(bahz)
library(rstan)
library(rethinking)
library(stringr)
options(mc.cores = parallel::detectCores())
source("src/functions.R")
library(tidyverse)
library(dplyr)
library(loo)
# Getting the minimal pre-compiled model to run -------------------
# Generate a dataset
data <- sim_geno_cline(transect_distances = seq(-300,300,20), n_ind = 30, Fis = 0,
                    decrease = T, center = 10, width = 80, pmin = 0.03, pmax = .95)
plot(x = data$transectDist, y = data$emp.p)
lines(x = data$transectDist, y = data$cline.p)
data2 <- rbind(data[1,])

make_prior_config(path = "~/Desktop/", name = "test.yaml")

prep_init_list("prior_config_template.yaml", tails = "left", chains = as.integer(3))

# Fit the model
fit_none <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                 type = "bi", tails = "none", chains = 3)
fit_left <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                      type = "bi", tails = "left", chains = 3)
fit_right <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                           type = "bi", tails = "right", chains = 3)
fit_mirror <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                           type = "bi", tails = "mirror", chains = 3)
fit_ind <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                            type = "bi", tails = "ind", chains = 3)

# ?nlist to get lists of prior stuff
cline_summary(fit_none)
cline_summary(fit_left)
cline_summary(fit_right)
cline_summary(fit_mirror)
cline_summary(fit_ind)


z1 <- loo::loo(fit_none, r_eff = relative_eff(fit_none))
z2 <- loo::loo(fit_left, r_eff = relative_eff(fit_left))
z3 <- loo::loo(fit_right, r_eff = relative_eff(fit_right))
z4 <- loo::loo(fit_mirror, r_eff = relative_eff(fit_mirror))
z5 <- loo::loo(fit_ind, r_eff = relative_eff(fit_ind))

loo::compare(z1, z2, z3, z4, z5)

# testing the new width check in prep_init_list
i <- 8000
z <- prep_init_list("prior_config_template.yaml", tails = "right", chains = as.integer(i))
widths <- rep(as.numeric(NA), times = i)
for (width in 1:i) {
  widths[width] <- z[[width]]$width
}
unique(unlist(widths) < 0) # Never makes a value below 0.


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


