
# Load packages -----------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(rstan)
#library(bahz)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayesplot)
library(loo)

# Generate phenotypic data from a cline -----------------------------------
x <-  seq(-200, 200, 20)

set.seed(123)
pheno <- sim_pheno_cline(transect_distances = x, n_ind = 20,
                         sigma = 6, decrease = T, center = 15, width = 30, pmin = 8, pmax = 22)
# Out of interest, simulate non-constant variance and see how it does
pheno <- sim_pheno_cline(transect_distances = x, n_ind = 20,
                         sigma = abs(rnorm(n = length(x), mean = 10, sd = 5)),
                         decrease = F, center = 150, width = 30, pmin = 8, pmax = 22)
# that worked fine. What about a more possibly interesting case: higher variance in the center
z <- ifelse(abs(15-x) <= 25, 14, 8)
pheno <- sim_pheno_cline(transect_distances = x, n_ind = 17,
                         sigma = z,
                         decrease = T, center = 15, width = 30, pmin = 35, pmax = 65)

site.means <- pheno %>%
  group_by(transectDist) %>%
  summarize(mean.pheno = mean(traitValue))
plot(pheno$transectDist, pheno$traitValue)
lines(x, site.means$mean.pheno, col = "red")
# So far, very simple: cline describes the mean phenotype, and variance is constant across the cline.



# Fit the model -----------------------------------------------------------

# Works! but no generated quantities yet.
z <- fit_pheno_cline(data = pheno,
                prior_file = "prior_config_template.yaml",
                chains = 4)

cline_summary(z)

plot_pheno_cline(z, data = pheno)
plot_pheno_cline(z, data = pheno, add.obs.pheno = T, point.col = "red", col = "blue")

pred <- predict_geno_cline(z, -200:200)



x <- cline_summary(z, show.all = T)
x$param

yrep <- as.data.frame(z) %>%
 select(starts_with("y_rep")) %>%
  as.matrix(.)

loo(z)


ll_mat <- as.data.frame(z) %>%
  select(starts_with("log_lik")) %>%
  as.matrix(.)

dim(yrep)

waic(ll_mat)
lm()


ppc_dens_overlay(y = pheno$traitValue, yrep = yrep[1:200,])

pred <- predict_geno_cline(z, distance = -200:200)

plot(pheno$transectDist, pheno$traitValue)
lines(pred$transectDist, pred$p, col = "blue")
