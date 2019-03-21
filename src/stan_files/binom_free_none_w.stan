#include /pre/license.stan

// Stan model to fit a sigmoid cline of allele frequencies,
// with no introgression tails.
// The minimum and maximum alele frequencies at the tails of the cline
// are estimated from the data

data {

#include /data/data_binom.stan

#include /priors/priors_binom.stan

}

parameters{
  real center; // the center of the cline, in km.
  real<lower=0> w; // the width of the cline. Also can't be negative
  real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
  real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
}

transformed parameters {
  vector[N] p; // a vector of expected p values for each site
    for ( i in 1:N ) { // for each site
    // calculate the predicted p
    if (decrease == 0) {
    p[i] =   pmin + (pmax - pmin) * (exp(w*(transectDist[i] - center))/(1 + exp(w * (transectDist[i] - center))));
    }
    if (decrease == 1) {
      p[i] =   pmin + (pmax - pmin) * (1-(exp(w*(transectDist[i] - center))/(1 + exp(w * (transectDist[i] - center)))));
    }
    }
}

model{
  // The statistical model
  pmax ~ uniform(p_l_max , p_u_max); // prior for pmax
  pmin ~ uniform(p_l_min, p_u_min); // prior for pmin
  w ~ exponential(p_scale_width); // prior for width.
  center ~ normal(p_center_1, p_center_2); //prior for center
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
