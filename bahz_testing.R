
# Load packages -----------------------------------------------------------
rm(list = ls())
library(devtools)
detach("package:bahz", unload=TRUE)
install_github("tjthurman/BAHZ", ref = "develop", auth_token = "c037ee28d7031c8d2cfb91dbefc686409bf05e68")
install.packages("~/Documents/Work/PhD:McGill/Projects/BAHZ/", repos = NULL, type = "source")
library(bahz)
library(rstan)
library(rethinking)
library(stringr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("src/functions.R")
library(tidyverse)

# Test out the pipeline ---------------------------------------------------
# Generate a dataset
data <- sim_data_from_cline(transect_distances = seq(0,600,20), n_ind = 40, Fis = 0.8,
                    decrease = F, center = 200, width = 38, pmin = 0.03, pmax = 1)
data2 <- data %>%
  mutate(transect_actual = transectDist) %>%
  mutate(transectDist = transectDist-mean(transectDist))
plot(data2$transectDist, data2$cline.p)
plot(data$transectDist, data$cline.p)

plot(data$transectDist, data$cline.p, type = "l", ylm = c(0,1))
lines(data$transectDist, data$emp.p, type = "l")

# Make a prior config file
make_prior_config(overwrite = F)

# Fit the model
bi_none <- fit_cline(data = data,
          prior_file = "prior_config_template.yaml",
          type = c("bi"),
          tails = c("none"),
          direction = c("inc"),
          chains = 3)

multi_none <- fit_cline(data = data,
                     prior_file = "prior_config_template.yaml",
                     type = "multi",
                     tails = "none",
                     direction = "inc")

precis(bi_none, prob = 0.05)
cline_summary(bi_none, prob = 0.05, method = "HPDI", show.all = F)


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


library(shinystan)
launch_shinystan(bi_left_limit)

bi_left@

precis(bi_left)
bi_right <- stan(
  model_code = models$bi_right_inc,
  data = data.list,
  chains = 3,
  iter = 10000,
  warmup = 3000,
  init = right.init.list
)
bi_mirror <- stan(
  model_code = models$bi_mirror_inc,
  data = data.list,
  chains = 3,
  iter = 10000,
  warmup = 3000,
  init = mirror.init.list
)
bi_ind <- stan(
  model_code = models$bi_ind_inc,
  data = data.list,
  chains = 3,
  iter = 10000,
  warmup = 3000,
  init = ind.init.list
)


bi_none@inits

plot(as.matrix(bi_none)[1:6000,1],type='l')




precis(bi_none)
cline_summary(bi_none, prob = .89, method = "HPDI")
precis(multi_none)
cline_summary(multi_none, prob = .97, method = "HPDI")

precis(bi_left)
cline_summary(bi_left, prob = .89)
compare(bi_none, bi_left)

precis(bi_right)
cline_summary(bi_right, prob = .89)

precis(bi_mirror)
cline_summary(bi_mirror, prob = .89)

precis(bi_ind)
cline_summary(bi_ind, prob = .89)

compare(bi_none, bi_left, bi_right, bi_mirror, bi_ind)



# library(rstan)
# library(rethinking)
# library(bahz)
# options(mc.cores = parallel::detectCores())
# # #
# set.seed(10)
# cline <- sim_data_from_cline(transect_distances = seq(0,600,20), n_ind = 40, Fis = 0,
#                             decrease = F, center = 300, width = 80)
# save(cline, file = "tests/testthat/reference_cline.Rdata")
# saveRDS(data, file = "tests/testthat/reference_cline")
# # #

# # Testing the centering
# # One strategy: take in the data however the user want to provide it,
# # center it before going to stan
# # and then modify the stan file after to put things back on the user scale
# mean.site <- mean(data$transectDist)
# data.c <- data
# data.c$transectDist <- data.c$transectDist - mean.site
#
# init_list <- list(list(width = .02),
#                    list(width = .02),
#                    list(width = .02),
#                    list(width = .02))
#
# stan_data <- load_cline_data(data.c, type = "bi")
#
# z <- stan(file = "src/stan_files/test_centered_c.stan", data = stan_data, init = init_list)
# z <- stan(file = "src/stan_files/test_centered_w.stan", data = stan_data)
# z <- stan(file = "src/stan_files/minimal.stan", data = stan_data, init = init_list)
#
# # problem seems to be that if the cline is too narrow, it has lots of trouble initializing
#
# z@inits
#
# cline_summary(z)
#
# 4/z@inits[[1]]$width
# 4/z@inits[[2]]$width
# 4/z@inits[[3]]$width
# 4/z@inits[[4]]$width
#
# summary(z)
#
# z@sim$samples[[1]][1]$center <- z@sim$samples[[1]][1]$center + mean.site
# z@sim$samples[[2]][1]$center <- z@sim$samples[[2]][1]$center + mean.site
# z@sim$samples[[3]][1]$center <- z@sim$samples[[3]][1]$center + mean.site
# z@sim$samples[[4]][1]$center <- z@sim$samples[[4]][1]$center + mean.site
#
# z@sim$samples[[1]][2]$width <- 4/z@sim$samples[[1]][2]$width
# z@sim$samples[[2]][2]$width <- 4/z@sim$samples[[2]][2]$width
# z@sim$samples[[3]][2]$width <- 4/z@sim$samples[[3]][2]$width
# z@sim$samples[[4]][2]$width <- 4/z@sim$samples[[4]][2]$width
#
#
# cline_summary(z)
# # this works, but only using precis, which re-calculates the stuff
# # The summary of the stanfit object is still the same, so I'd have to recalculate it somehow.
# # The problem is when initiation values are way too small: just use half the distance of the transect?
#
#
# cline_summary(z)
#
#
# # Option 2: do it with transformed parameters and/or generated quantities somehow.
#
# stan_data2 <-  with(data, list(N = length(transectDist),
#                                nFocalAllele = 2*AA + Aa,
#                                nTotalAlleles = 2*N,
#                                transectDist = transectDist,
#                                meanDist = mean(transectDist)))
# x <- stan(file = "src/stan_files/test_centered_c_instan.stan", data = stan_data2, init = init_list)
# precis(x)
# cline_summary(x)
#
# # This works, but has both c and center in the final stanfit object, which is a bit annoying.
#
#
# # After some fussing things, seems like the first method is easier. Can just build it in to the main
# # model fitting wrapper
#
#
#
# # testing Fitting with pre-compiled code
# init_list <- make_init_list("inst/extdata/prior_config_template.yaml", tails = "none", chains = as.integer(3))
# z <- bahz::test_fit_cline(stan_data, init_list)
#
# ws <- rexp(100, 10)
#
# plot(density(4/ws))
