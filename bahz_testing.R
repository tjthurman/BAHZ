
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
                    decrease = F, center = 0, width = 20, pmin = 0.03, pmax = .95)
# Import to stan using function I already made:
st_dat <- load_cline_data(data, type = "bi")
make_prior_config()
# Make a list of start values:
init <- make_init_list("prior_config_template.yaml", tails = "none", chains = as.integer(3))

init2 <- init
init2[[1]]$width <- 4/20
init2[[2]]$width <- 4/20.4
init2[[3]]$width <- 4/21
names(init2[[1]])[2] <- "w"
names(init2[[2]])[2] <- "w"
names(init2[[3]])[2] <- "w"

priors <- list(p_m_center = 0,
               p_sd_center = 150,
               p_m_width = 20,
               p_sd_width = 100,
               p_l_min = 0,
               p_u_min = 0.2,
               p_l_max = 0.8,
               p_u_max = 1,
               p_scale_width = 10)
priors2 <- list(p_m_center = 0,
               p_sd_center = 150,
               p_scale_width = 10,
               p_l_min = 0,
               p_u_min = 0.2,
               p_l_max = 0.8,
               p_u_max = 1,
               p_m_width =0,
               p_sd_width = 0)

z_p <- test_fit_cline(stan_data = c(st_dat, priors), init_list = init, chains = 3, model = "binom_free_none_all_priors")

z_p2 <- test_fit_cline(stan_data = c(st_dat, priors2), init_list = init2, chains = 3, model = "binom_free_none_all_priors_w")


# ?nlist to get lists of prior stuff
cline_summary(z_p)
cline_summary(z_p2)

mean(4/as.data.frame(z_p2)$w)

# For diagnoses when things go wrong, but takes forever
# pdf("pairs.pdf", width = 20, height = 20)
# pairs.default(z_p2)
# dev.off()


# Testing centering of center --------------------------------------------
none.init.list.actual <- list(list(center = 285, width = 38, pmin = 0.05, pmax = .99))
bi_none_actual <- stan(
  file = "src/stan_models/binomial/binom_free_none_nolim.stan",
  data = load_cline_data(data, type = "bi"),
  chains = 1,
  iter = 10000,
  warmup = 3000,
  init = none.init.list.actual
)
precis(bi_none_actual)

none.init.list.center <- list(list(width = 38))
bi_none_center <- stan(
  file = "src/stan_models/binomial/binom_free_none_nolim.stan",
  data = load_cline_data(data2, type = "bi"),
  chains = 1,
  iter = 10000,
  warmup = 3000,
  init = none.init.list.center
)
precis(bi_none_center)
bi_none_center@sim$samples[[1]]$center <- bi_none_center@sim$samples[[1]]$center + mean(data2$transect_actual)
# Centering by subtracting the mean from the collecting sites seems to work well!
# Can add the mean back on after estimation to put it back on user scale.

# Can also use the centering thing on priors: the sd will stay constant.
z <- rnorm(1000, mean = 200, sd = 2)
center.z <- z - mean(z)
mean(center.z)
sd(z)
sd(center.z)


# Will stan run with no priors? -------------------------------------------
bi_none_nopriors <- stan(
  file = "src/stan_models/binomial/binom_free_none_nolim_nopriors.stan",
  data = load_cline_data(data, type = "bi"),
  chains = 3,
  iter = 10000,
  warmup = 3000,
  init = none.init.list
)

precis(bi_none_nopriors)

none.init.list <- list(list(center = 300, width = 38, pmin = 0.05, pmax = .99),
list(center = 357, width = 28, pmin = 0.15, pmax = .87),
list(center = 290, width = 18, pmin = 0.02, pmax = .8))

# Yes, stan will run with no priors! I think the default, then, is that it
# puts uniform -inf, Inf priors, which against the whole point of using stan,
# really. Still, good to test.


# Re-scaling of width? ----------------------------------------------------

# Is it possible to rescale the width parameter as well, using 4/width

# Testing on how to convert from normal distribution on width scale
# to normal distribution on w scale.
m <- 100
sd <- 30
CV <- sd/m

new.mean <- 4/m
new.sd <- new.mean*CV

new.sd/new.mean == CV

width <- rnorm(1000, mean = m, sd = sd)
w <- rnorm(1000, mean = new.mean, sd= new.sd)

w2 <- 4/width

w.to.width <- 4/w
range(width)
range(w.to.width)
mean(width)
mean(w.to.width)
mean(w2)
sd(width)
sd(w.to.width)
sd(w2)
dens(4/width)
dens(w2)

# After a lot of testing, the answer seems to be no.
# It would be nice to re-scale width to ease the starting and allow it
# to start from default values. But, I'm nt sure there's any way to convert
# the user-supplied prior to a sensible prior on the 4/width scale.

# What parameters really need to be started? ------------------------------






# Does the delta limit do anything? ---------------------------------------
left.init.list <- list(list(center = 285, width = 38),
                       list(center = 285, width = 38),
                       list(center = 285, width = 38))
bi_none_limit <- stan(
  file = "src/stan_models/binomial/binom_free_none.stan",
  data =  load_cline_data(data, type = "bi"),
  chains = 3,
  iter = 10000,
  warmup = 3000,
  init = left.init.list, verbose = T
)

bi_left_limit <- stan(
  file = "src/stan_models/binomial/binom_free_left.stan",
  data =  load_cline_data(data, type = "bi"),
  chains = 3,
  iter = 10000,
  warmup = 3000,
  init = left.init.list, verbose = T
)

bi_left_nolimit <- stan(
  file = "src/stan_models/binomial/binom_free_left_nodeltalim.stan",
  data =  load_cline_data(data, type = "bi"),
  chains = 3,
  iter = 10000,
  warmup = 3000,
  init = left.init.list, verbose = T
)
# When simulating data with a deltaL of 8, didn't matter much.
# Let's do it with no tail.
cline_summary(bi_left_limit)
cline_summary(bi_left_nolimit)

# When simulating data wihtout a tail, the limit made a difference,
# keeping the estimate of delta much further down (by 50%)
compare(bi_left_limit, bi_left_nolimit)
compare(bi_left_nolimit, bi_none_limit)
# The WAIC is better for bi_none, it takes about 64% of the weight.

