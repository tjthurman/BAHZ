data{
  int<lower=1> N; // number of sites sampled
  int nFocalAllele[N]; // number of alleles of the focal type, in this case hydara alleles
  int nTotalAlleles[N]; // total number of alleles sampled (2*number of diploid individuals)
  real transectDist[N]; // distance along transect
}

parameters{
  real center;
  real<lower=0> width; // the width of the cline. Also can't be negative
}

transformed parameters {
  vector[N] p; // a vector of expected p values for each site
    for ( i in 1:N ) { // for each site
    // calculate the predicted p
    p[i] = (exp(width*(transectDist[i] - center))/(1 + exp(width * (transectDist[i] - center))));
  }
}
model{
  // The statistical model
  width ~ exponential(1); // prior for width.
  center ~ normal(0, 10); //prior for center
  // and the likelihood: observed allele counts follow a binomial liklihood,
  // based on the number of alleles sampled and the estimated allele frequency.
  nFocalAllele ~ binomial(nTotalAlleles, p);
}
// generated quantities{
//   // these are used mostly for calculating WAIC and DIC
//   vector[N] p; // calculate a vector of expected p values
//   real dev; // calculate a deviance value for the model overall based on the draws from the posterior
//   vector[N] log_lik; // calculate a vector of log-liklihoods for each observed allele count based on the draws from the posterior
//    vector[N] y_rep; // for posterior preictive checks, an allele count drawn from a binomial distribution
//   for ( i in 1:N ) { // for each site
//     // calculate the predicted p for each site.
//     p[i] = pmin + (pmax - pmin) * (exp(width*(transectDist[i] - center))/(1 + exp(width * (transectDist[i] - center))));
//     // and calculate the log likelihood of the observed allele counts, given the predicted p
//     log_lik[i] = binomial_lpmf(nFocalAllele[i] | nTotalAlleles[i], p[i]);
//     // and calculate an expected number of focal alleles, given the same number of samples
//     y_rep[i] = binomial_rng(nTotalAlleles[i], p[i]);
//   }
//   // calculate deviance of the model
//   dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);
// }


