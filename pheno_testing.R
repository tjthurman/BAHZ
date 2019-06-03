
# Load packages -----------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(rstan)
#library(bahz)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Generate phenotypic data from a cline -----------------------------------
x <-  seq(0, 300, 20)

set.seed(123)
sim_geno_cline(10, 10, 0, T, 10, 10)
pheno <- sim_pheno_cline(transect_distances = x, n_ind = 20,
                         sigma = 6, decrease = F, center = 150, width = 30, pmin = 8, pmax = 22)
# Out of interest, simulate non-constant variance and see how it does
pheno <- sim_pheno_cline(transect_distances = x, n_ind = 20,
                         sigma = abs(rnorm(n = length(x), mean = 10, sd = 5)),
                         decrease = F, center = 150, width = 30, pmin = 8, pmax = 22)
# that worked fine. What about a more possibly interesting case: higher variance in the center
z <- ifelse(abs(150-x) <= 15, 14, 7)
pheno <- sim_pheno_cline(transect_distances = x, n_ind = 17,
                         sigma = z,
                         decrease = F, center = 150, width = 30, pmin = 35, pmax = 65)

site.means <- pheno %>%
  group_by(transectDist) %>%
  summarize(mean.pheno = mean(traitValue))
plot(pheno$transectDist, pheno$traitValue)
lines(x, site.means$mean.pheno, col = "red")
# So far, very simple: cline describes the mean phenotype, and variance is constant across the cline.

# Can I fix the bad vectorization -----------------------------------------

# Trying a model based on the Ragged array design in
# section 8.2 of the Stan User's guide v2.19
prelim_pheno_stan2 <- "
data{
int<lower=0> N; // number of individuals
int<lower=0> K; // number of sites
vector[N] pheno; // phenotype for an individual
int s[K]; //number of individuals sampled per site
real transectDist[K]; // distance along transect for a site


}
parameters{
real center; // the center of the cline, in km.
real<lower=0> width; // the width of the cline. Also can't be negative
real pmin; // minimum mean pheno value
real pmax; // maximum mean pheno value
real<lower =0> sigma; // phenotypic variance (constant for now, must be positive)
}

transformed parameters{
vector[K] p; // the expected phenotype for each site
for (i in 1:K)
{
  p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
}
}

model{
int pos;
pos = 1;
sigma ~ normal(0, 100);
pmin ~ normal(50, 100);
pmax ~ normal(50, 100);
width ~ normal(50, 100);
center ~ normal(100, 100);

for (k in 1:K) {
segment(pheno, pos, s[k]) ~ normal(p[k], sigma);
pos = pos + s[k];
}
}

"


prep_pheno_data(dplyr::sample_frac(pheno, 1))

z2 <- stan(model_code = prelim_pheno_stan2,
           data = prep_pheno_data(dplyr::sample_frac(pheno, 1)),
           init = list(list(center = rnorm(1,100, 100),
                            width = abs(rnorm(1, 50, 100)),
                            pmin = abs(rnorm(1, 50, 100)),
                            pmax = abs(rnorm(1, 50,100)),
                            sigma = abs(rnorm(1, 0,100))),
                       list(center = rnorm(1,100, 100),
                            width = abs(rnorm(1, 50, 100)),
                            pmin = abs(rnorm(1, 50, 100)),
                            pmax = abs(rnorm(1, 50,100)),
                            sigma = abs(rnorm(1, 0,100))),
                       list(center = rnorm(1,100, 100),
                            width = abs(rnorm(1, 50, 100)),
                            pmin = abs(rnorm(1, 50, 100)),
                            pmax = abs(rnorm(1, 50,100)),
                            sigma = abs(rnorm(1, 0,100))),
                       list(center = rnorm(1,100, 100),
                            width = abs(rnorm(1, 50, 100)),
                            pmin = abs(rnorm(1, 50, 100)),
                            pmax = abs(rnorm(1, 50,100)),
                            sigma = abs(rnorm(1, 0,100)))))
cline_summary(z)
cline_summary(z2)

pred2 <- predict_geno_cline(z2, distance = 0:200)

plot(inds$distance, inds$pheno)
lines(pred$transectDist, pred$p, col = "red")
lines(pred2$transectDist, pred2$p, col = "blue")

# That worked! same numerical result, and the stan fit object is way smaller.
cline_summary(z2, show.all = T)
# Only found mean phenotypes for the 21 sites, not for all
# 201 individuals
cline_summary(z, show.all = T)
