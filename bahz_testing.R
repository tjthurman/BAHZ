
# Load packages -----------------------------------------------------------
rm(list = ls())
library(devtools)
detach("package:bahz", unload=TRUE)
install_github("tjthurman/BAHZ", ref = "feature/rstantools_attempt", auth_token = "c037ee28d7031c8d2cfb91dbefc686409bf05e68")
install.packages("~/Documents/Work/PhD:McGill/Projects/BAHZ/", repos = NULL, type = "source")
library(bahz)
library(rstan)
library(rethinking)
library(stringr)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("src/functions.R")
library(tidyverse)
library(dplyr)
# Getting the minimal pre-compiled model to run -------------------
# Generate a dataset
data <- sim_data_from_cline(transect_distances = seq(-300,300,20), n_ind = 40, Fis = 0,
                    decrease = T, center = 0, width = 50, pmin = 0.03, pmax = .95)

data2 <- rbind(data[1,])

plot(data$transectDist, data$emp.p)
# Import to stan using function I already made:
st_dat <- load_cline_data(data, type = "bi")
st_dat <- load_cline_data(data, type = "multi")
make_prior_config()
# Make a list of start values:
init <- make_init_list("prior_config_template.yaml", tails = "none", chains = as.integer(3))
prior_list <- make_prior_list("prior_config_template.yaml")

init2 <- init
init2[[1]]$width <- 4/20
init2[[2]]$width <- 4/20.4
init2[[3]]$width <- 4/21
names(init2[[1]])[2] <- "w"
names(init2[[2]])[2] <- "w"
names(init2[[3]])[2] <- "w"

priors <- list(p_center_1 = 0,
               p_center_2 = 150,
               p_dist_center = 0,
               p_m_width = 20,
               p_sd_width = 100,
               p_l_min = 0,
               p_u_min = 0.2,
               p_l_max = 0.8,
               p_u_max = 1,
               p_scale_width = 0)
priors2 <- list(p_m_center = 0,
               p_sd_center = 150,
               p_scale_width = 10,
               p_l_min = 0,
               p_u_min = 0.2,
               p_l_max = 0.8,
               p_u_max = 1,
               p_m_width =0,
               p_sd_width = 0)

z_p <- test_fit_cline(stan_data = c(st_dat, priors), init_list = init, chains = 3, model = "binom_free_none_width")

z_p2 <- test_fit_cline(stan_data = c(st_dat, priors2), init_list = init2, chains = 3, model = "binom_free_none_w")


# ?nlist to get lists of prior stuff
cline_summary(z_p)
cline_summary(z_p2)

mean(4/as.data.frame(z_p2)$w)


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


