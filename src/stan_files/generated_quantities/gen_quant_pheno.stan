// these are used mostly for calculating WAIC and DIC
// real dev; // calculate a deviance value for the model overall based on the draws from the posterior
//real log_lik[N*k]; // calculate a vector of log-liklihoods for each observed phenotype based on the draws from the posterior
// log_lik = rep_array(0, N);
real y_rep[N]; // for posterior preictive checks, a phenotype drawn from the posterior
for (k in 1:K) { //for each site
  //log_lik[starts[k]:ends[k]] = normal_lpdf(pheno[starts[k]:ends[k]] | p[k], sigma);
  y_rep[starts[k]:ends[k]] = normal_rng(rep_array(p[k], s[k]), rep_array(sigma, s[k]));
  }
// calculate deviance of the model
// NEED TO FIGURE OUT HOW TO CALCULATE FOR THE PHENO STUFF
// dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);

