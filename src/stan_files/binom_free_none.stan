#include /pre/license.stan

// Stan model to fit a sigmoid cline of allele frequencies,
// with no introgression tails.
// The minimum and maximum alele frequencies at the tails of the cline
// are estimated from the data

#include /data_binom.stan

parameters{
  real<lower=0> center; // the center of the cline, in km. Can't be negative
  real<lower=0> width; // the width of the cline. Also can't be negative
  real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
  real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
}
model{
  // The statistical model
  vector[N] p; // a vector of expected p values for each site
  pmax ~ uniform(0.8 , 1); // prior for pmax
  pmin ~ uniform(0 , 0.2); // prior for pmin
  width ~ normal(50, 100); // prior for width.
  center ~ normal(350, 150); //prior for center
  for ( i in 1:N ) { // for each site
    // calculate the predicted p
    p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
  }
  // and the likelihood: observed allele counts follow a binomial liklihood,
  // based on the number of alleles sampled and the estimated allele frequency.
  nFocalAllele ~ binomial(nTotalAlleles, p);
}
generated quantities{
  // these are used mostly for calculating WAIC and DIC
  vector[N] p; // calculate a vector of expected p values
  real dev; // calculate a deviance value for the model overall based on the draws from the posterior
  vector[N] log_lik; // calculate a vector of log-liklihoods for each observed allele count based on the draws from the posterior
   vector[N] y_rep; // for posterior preictive checks, an allele count drawn from a binomial distribution
  for ( i in 1:N ) { // for each site
    // calculate the predicted p for each site.
    p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
    // and calculate the log likelihood of the observed allele counts, given the predicted p
    log_lik[i] = binomial_lpmf(nFocalAllele[i] | nTotalAlleles[i], p[i]);
    // and calculate an expected number of focal alleles, given the same number of samples
    y_rep[i] = binomial_rng(nTotalAlleles[i], p[i]);
  }
  // calculate deviance of the model
  dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);
}
