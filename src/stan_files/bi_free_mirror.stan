#include /pre/license.stan

data {

#include /data/data_binom.stan
#include /priors/priors_all.stan
#include /priors/priors_mirror_tail.stan


}

parameters{
  real center; // the center of the cline, in km.
  real<lower=0> width; // the width of the cline. Also can't be negative
  real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
  real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
  real<lower=0> deltaM; // the location of the tails.
  real<lower=0,upper=1> tauM; // the ratio of the slope of the exponential tail to the slope of the sigmoid center. Must be between 0 and 1.
}

transformed parameters {
  vector[N] p; // a vector of expected p values for each site
  for ( i in 1:N ) { // for each site
    // calculate the predicted p
     if (transectDist[i] <= center - deltaM) { // if the site is in the left tail
        // predicted value of p is determined by the exponential equation for the left tail
        if (decrease == 0) {
           p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaM/width)))*exp((4*tauM*(transectDist[i] - center + deltaM)/width)/(1 + exp(-4*deltaM/width)));
        }
        if (decrease == 1) {
           p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaM/width)))*exp((4*tauM*(transectDist[i] - center + deltaM)/width)/(1 + exp(-4*deltaM/width))));
        }
     } else if (transectDist[i] >= center + deltaM) { // if the site is in the right tail
        // predicted value of p is determined by the exponential equation for the right tail
        if (decrease == 0) {
          p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaM/width)))*exp((-4*tauM*(transectDist[i] - center - deltaM)/width)/(1 + exp(-4*deltaM/width))));
        }
        if (decrease == 1) {
          p[i] =  pmin + (pmax-pmin)*((1/(1 + exp(4*deltaM/width)))*exp((-4*tauM*(transectDist[i] - center - deltaM)/width)/(1 + exp(-4*deltaM/width))));
        }
     } else { // otherwise, the site is in the exponential center
    // so use the regular sigmoid cline equation
          if (decrease == 0) {
              p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
          }
          if (decrease == 1) {
              p[i] =   pmin + (pmax - pmin) * (1-(exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width))));
          }
        }
    }
}

model{
  // The statistical model

#include /model/model_mirror_tail.stan
#include /model/model_ps.stan
#include /model/model_width.stan
#include /model/model_center.stan

  // and the likelihood: observed allele counts follow a binomial liklihood,
  // based on the number of alleles sampled and the estimated allele frequency.
  nFocalAllele ~ binomial(nTotalAlleles, p);
}
generated quantities{
  // these are used mostly for calculating WAIC and DIC
  real dev; // calculate a deviance value for the model overall based on the draws from the posterior
  vector[N] log_lik; // calculate a vector of log-liklihoods for each observed allele count based on the draws from the posterior
  vector[N] y_rep; // for posterior preictive checks, an allele count drawn from a binomial distribution
  for ( i in 1:N ) { // for each site
    // and calculate the log likelihood of the observed allele counts, given the predicted p
    log_lik[i] = binomial_lpmf(nFocalAllele[i] | nTotalAlleles[i], p[i]);
    // and calculate an expected number of focal alleles, given the same number of samples
    y_rep[i] = binomial_rng(nTotalAlleles[i], p[i]);
  }
  // calculate deviance of the model
  dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);
}
