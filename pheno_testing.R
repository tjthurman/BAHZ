
# Load packages -----------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(rstan)
#library(bahz)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# Generate phenotypic data from a cline -----------------------------------
x <-  seq(0, 200, 30)

set.seed(123)
cline <- general_cline_eqn(transectDist = x,
                  decrease = F,
                  center = 100,
                  width = 30, pmin = 10, pmax = 40)
inds <- NULL
for (site in 1:length(x)) {
  site.res <- data.frame(site = site,
                         distance = x[site],
                         pheno = rnorm(10, mean = cline[site], sd = 5))
  inds <- rbind(inds, site.res)
}

plot(inds$distance, inds$pheno)
lines(x, cline)
means <- inds %>%
  group_by(distance) %>%
  summarize(m = mean(pheno))
lines(means$distance, means$m, col = "red")
# So far, very simple: cline describes the mean phenotype, and variance is constant across the cline.

# A simple cline model
prelim_pheno_stan <- "
data{
int N; // number of individuals sampled
vector[N] pheno; // phenotype for an individual
vector[N] transectDist; // distance along transect
}
parameters{
real center; // the center of the cline, in km.
real<lower=0> width; // the width of the cline. Also can't be negative
real pmin; // minimum mean pheno value
real pmax; // maximum mean pheno value
real<lower =0> sigma; // phenotypic variance (constant for now, must be positive)
}

transformed parameters{
vector[N] p; // the expected phenotype for each individual
for (i in 1:N)
{
p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
}
}

model{
  sigma ~ normal(0, 100);
  pmin ~ normal(50, 100);
  pmax ~ normal(50, 100);
  width ~ normal(50, 100);
  center ~ normal(100, 100);
  pheno ~ normal(p, sigma);
}

"
# Does it work?
z <- stan(model_code = prelim_pheno_stan, data = list(N = dim(inds)[1],
                                                 pheno = inds$pheno,
                                                 transectDist = inds$distance))
# yes!
cline_summary(z)
# Only problem: not vectorizing well. Its calculating and sampling a separate p for every data point,
# which we do not want it to do. We'd rather it calculate a p for every site. Might have to structure the data
# differently

pred <- predict_geno_cline(z, distance = 0:200)

plot(inds$distance, inds$pheno)
lines(pred$transectDist, pred$p, col = "red")

