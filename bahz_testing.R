
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



# Testing geno models -----------------------------------------------------
# Generate a dataset
set.seed(22)
data <- sim_geno_cline(transect_distances = seq(-300,300,20), n_ind = 40, Fis = .9,
                       decrease = F, center = 10, width = 35, pmin = 0.08, pmax = .95, deltaR = 12, tauR = 0.25)

plot(x = data$transectDist, y = data$emp.p)
lines(x = data$transectDist, y = data$cline.p)



none_bi <- fit_geno_cline(data = data2, prior_file = "~/Desktop/geno_priors.yaml",
                              type = "bi", tails = "none", ignore_data = F)
none_multi <- fit_geno_cline(data = data2, prior_file = "~/Desktop/geno_priors.yaml",
                              type = "multi", tails = "none", ignore_data = F)

left_bi <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                               type = "bi", tails = "left", ignore_data = F)
left_multi <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                           type = "multi", tails = "left", ignore_data = F)

right_bi <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                           type = "bi", tails = "right", ignore_data = F)
right_multi <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                            type = "multi", tails = "right", ignore_data = F)

mirror_bi <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                            type = "bi", tails = "mirror", ignore_data = F)
mirror_multi <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                             type = "multi", tails = "mirror", ignore_data = F)

ind_bi <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                             type = "bi", tails = "ind", ignore_data = F)
ind_multi <- fit_geno_cline(data = data, prior_file = "~/Desktop/geno_priors.yaml",
                          type = "multi", tails = "ind", ignore_data = F)



rethinking::compare(none_bi, left_bi, right_bi, mirror_bi, ind_bi)
cline_summary(none_bi)
cline_summary(none_bi, show.all = T)
cline_summary(none_multi)
cline_summary(left_bi)
cline_summary(left_multi)
cline_summary(right_bi)
cline_summary(right_multi)
cline_summary(mirror_bi)
cline_summary(mirror_multi)
cline_summary(ind_bi)
cline_summary(ind_multi)

# Prior predictive checks
library(bayesplot)
ppc


y_rep_bi <- as.matrix(none_bi, pars = names(none_bi)[which(str_detect(names(none_bi), pattern = "^y_rep"))])
y_rep_multi <- as.matrix(none_multi, pars = names(none_multi)[which(str_detect(names(none_multi), pattern = "^y_rep"))])


ppc_intervals(y = data$AA*2 + data$Aa, yrep = y_rep_bi, prob = 0.99)


# Testing pheno models ----------------------------------------------------

set.seed(839)
data <- sim_pheno_cline(transect_distances = seq(-200, 150, length.out = 12), n_ind = as.integer(rnorm(n = 12, mean = 30, sd = 7)),
                        sigma = abs(rnorm(n = 12, mean = 40, sd = 8)), decrease = T, center = 9, width = 72, pmin = 220, pmax = 762)

constant <- fit_pheno_cline(data, prior_file = "~/Desktop/pheno_priors.yaml", pheno_variance = "constant", ignore_data = F)
independent <- fit_pheno_cline(data, prior_file = "~/Desktop/pheno_priors.yaml", pheno_variance =  "independent", ignore_data = F)
pooled <- fit_pheno_cline(data, prior_file = "~/Desktop/pheno_priors.yaml", pheno_variance = "pooled", ignore_data = F)

names(constant)
cline_summary(constant, show.all = F)
cline_summary(constant, show.all = T)

cline_summary(independent, show.all = F)
cline_summary(pooled)



# Some notes on prior predictive checks -----------------------------------
# Note that for parameters with limits (like width), values
# that fall outside the allowed threshold are ignored.

# So, if the prior for center is normal with mean 50 and sd 50, the
# actual prior distribution in the end is:

z <- rnorm(n = 1000, mean = 50, sd = 50)
mean(z[z>0]) # 61, not
mean(z) # ~50

# To get bulk and tail ESS ------------------------------------------------
# A little odd and annoying that these functions don't just work on
# stanfit objects.
# But, easiest thing to do is use the monitor function in rstan.

# Can run the function on the stanfit object by itself,
# Which apparently isn't supposed to work.
monitor(none_bi)

# This sort of works, but gets the number of warmup iterations wrong also:
monitor(as.array(none_bi), warmup = 0)

# This method seems to work properly:
monitor(extract(none_bi, permuted = F, inc_warmup = T))


# Commented on a Github issue about this.

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


