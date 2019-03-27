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
