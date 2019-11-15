
# Load packages -----------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(rstan)
#library(bahz)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayesplot)
library(loo)
library(rstanarm)

# Generate phenotypic data from a cline -----------------------------------

x <- seq(-100, 100, 20)
set.seed(123)
n.sites <- 10
tenSites <- sim_pheno_cline(transect_distances = seq(from = -100, to = 100, length.out = n.sites),
                            n_ind = 100/n.sites,
                            sigma = 2, decrease = F, center = 15, width = 30, pmin = 8, pmax = 22)
n.sites <- 5
fiveSites <- sim_pheno_cline(transect_distances = seq(from = -100, to = 100, length.out = n.sites),
                                      n_ind = 100/n.sites,
                                      sigma = 2, decrease = F, center = 15, width = 30, pmin = 8, pmax = 22)
n.sites <- 20
twentySites <- sim_pheno_cline(transect_distances = seq(from = -100, to = 100, length.out = n.sites),
                            n_ind = 100/n.sites,
                            sigma = 2, decrease = F, center = 15, width = 30, pmin = 8, pmax = 22)

# Out of interest, simulate non-constant variance and see how it does
pheno <- sim_pheno_cline(transect_distances = x, n_ind = 20,
                         sigma = 400,
                         decrease = F, center = 15, width = 30, pmin = 1200, pmax = 2400)
pheno <- sim_pheno_cline(transect_distances = x, n_ind = as.integer(abs(rnorm(length(x), mean = 15, sd = 9))),
                         sigma = rnorm(length(x), mean = 2, sd = 0.75),
                         decrease = F, center = 15, width = 30, pmin = 12, pmax = 24)
# that worked fine. What about a more possibly interesting case: higher variance in the center
z <- ifelse(abs(15-x) <= 30, 15, 4)
z <- ifelse(abs(15-x) <= 30, 15, 4)
pheno <- sim_pheno_cline(transect_distances = x, n_ind = 17,
                         sigma = z,
                         decrease = T, center = 15, width = 30, pmin = 5, pmax = 28)

site.means <- pheno %>%
  group_by(transectDist) %>%
  summarize(mean.pheno = mean(traitValue))
plot(pheno$transectDist, pheno$traitValue)
lines(x, site.means$mean.pheno, col = "red")
# So far, very simple: cline describes the mean phenotype, and variance is constant across the cline.



# Fit the model -----------------------------------------------------------

# Works! but no generated quantities yet.
constant_sig <- fit_pheno_cline(data = pheno,
                prior_file = "prior_config_template.yaml", pheno_variance = "constant",
                chains = 4)
pool_sig <- fit_pheno_cline(data = pheno,
                            prior_file = "prior_config_template.yaml", pheno_variance = "pooled",
                            chains = 4)
ind_sig <- fit_pheno_cline(data = pheno,
                     prior_file = "prior_config_template.yaml", pheno_variance = "independent",
                     chains = 4)

cline_summary(constant_sig, show.all = F)
cline_summary(pool_sig, show.all = F)
cline_summary(ind_sig, show.all = F)

cline_summary(constant_sig, show.all = T)
cline_summary(pool_sig, show.all = T)
cline_summary(ind_sig, show.all = T)



plot_pheno_cline(ind_sig, data = pheno)
plot_pheno_cline(ind_sig, data = pheno, add.obs.pheno = T, point.col = "red", col = "blue")

plot_pheno_cline(constant_sig, data = pheno)
plot_pheno_cline(constant_sig, data = pheno, add.obs.pheno = T, point.col = "red", col = "blue")
plot_pheno_cline(pool_sig, data = pheno, add.obs.pheno = T, point.col = "red", col = "blue")


waic_constant <- waic(extract_log_lik(constant_sig))
waic_pool <- waic(extract_log_lik(pool_sig))
waic_ind <- waic(extract_log_lik(ind_sig))

compare(waic_pool, waic_constant,  waic_ind)

loo_constant <- loo(constant_sig)
loo_pool <- loo(pool_sig)
loo_ind <- loo(ind_sig)

compare(loo_pool, loo_constant,  loo_ind)

pred <- predict_geno_cline(z, -200:200)

pp <- as.data.frame(constant_sig) %>%
  select(starts_with("y_rep"))

ppc_hist(y = pheno$traitValue, yrep = as.matrix(pp))

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


five <- fit_pheno_cline(data = fiveSites,
                                prior_file = "prior_config_template.yaml", pheno_variance = "constant",
                                chains = 4)

ten <- fit_pheno_cline(data = tenSites,
                                prior_file = "prior_config_template.yaml", pheno_variance = "constant",
                                chains = 4)

twenty <- fit_pheno_cline(data = twentySites,
                       prior_file = "prior_config_template.yaml", pheno_variance = "constant",
                       chains = 4)
cline_summary(five)
cline_summary(ten)
cline_summary(twenty)
plot_pheno_cline(five, fiveSites, add.obs.pheno = T)
plot_pheno_cline(ten, tenSites, add.obs.pheno = T)
plot_pheno_cline(twenty, twentySites, add.obs.pheno = T)
