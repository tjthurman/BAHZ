// these are used mostly for calculating WAIC and DIC
//real dev; // calculate a deviance value for the model overall based on the draws from the posterior
vector[N] log_lik; // calculate a vector of log-liklihoods for each observed genotype count based on the draws from the posterior
 int y_rep[N,3]; // for posterior preictive checks, a genotype count drawn from a multinomial distribution
for ( i in 1:N ) { // for each site
  // and calculate the log likelihood of the observed allele counts, given the predicted p
  log_lik[i] = multinomial_lpmf(genos[i] | probs[i]);
  // and calculate an expected number of focal alleles, given the same number of samples
  y_rep[i] = multinomial_rng(probs[i], sum(genos[i]));
}
 // calculate deviance of the model
 // isn't working like this, stan doesn't like the specifications for geno and probs
 // need to figure out, or see if I can just do it post-hoc from the log likelihood
 //dev = (-2)*multinomial_lpmf(genos | probs);
