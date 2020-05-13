// these are used mostly for calculating WAIC and DIC
// real dev; // calculate a deviance value for the model overall based on the draws from the posterior
real log_lik[N]; // calculate a vector of log-liklihoods for each observed phenotype based on the draws from the posterior
real y_rep[N]; // for posterior preictive checks, a phenotype drawn from the posterior
for (k in 1:K) { //for each site
  y_rep[starts[k]:ends[k]] = normal_rng(rep_array(p[k], s[k]), rep_array(site_sigma[k], s[k]));
  for (n in 1:N) { //for each individual
    if (indDist[n] == transectDist[k]) {
      log_lik[n] = normal_lpdf(pheno[n] | p[k], site_sigma[k]);
    }
  }
}
// calculate deviance of the model
// NEED TO FIGURE OUT HOW TO CALCULATE FOR THE PHENO STUFF
// dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);

