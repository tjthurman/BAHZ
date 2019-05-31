
# Load packages -----------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(rstan)
#library(bahz)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# Generate phenotypic data from a cline -----------------------------------
x <-  seq(0, 300, 5)

set.seed(123)
cline <- general_cline_eqn(transectDist = x,
                  decrease = F,
                  center = 150,
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
                                                 transectDist = inds$distance),
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
# yes!
cline_summary(z, show.all = T)
# Only problem: not vectorizing well. Its calculating and sampling a separate p for every data point,
# which we do not want it to do. We'd rather it calculate a p for every site. Might have to structure the data
# differently

pred <- predict_geno_cline(z, distance = 0:200)

plot(inds$distance, inds$pheno)
lines(pred$transectDist, pred$p, col = "red")


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

n.per.site <- inds %>%
  group_by(distance) %>%
  tally()
z2 <- stan(model_code = prelim_pheno_stan2,
           data = list(N = dim(inds)[1],
                       K = length(unique(inds$distance)),
                       pheno = inds$pheno,
                       s = n.per.site$n,
                       transectDist = unique(inds$distance)),
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
