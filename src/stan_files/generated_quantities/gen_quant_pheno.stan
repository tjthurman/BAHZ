// these are used mostly for calculating WAIC and DIC
// real dev; // calculate a deviance value for the model overall based on the draws from the posterior
// vector[N] log_lik; // calculate a vector of log-liklihoods for each observed phenotype based on the draws from the posterior
// vector[N] y_rep; // for posterior preictive checks, an phenotype drawn from the posterior
// int poss;
// poss = 1;
//
// for (k in 1:K) { //for each site
//   segment(log_lik, poss, s[k]) = normal_lpdf(segment(pheno, poss, s[k]) | p[k], sigma);
//   segment(y_rep, pos2, s[k]) = normal_rng(p[k], sigma);
//   poss = poss + s[k];
//   }
// calculate deviance of the model
// NEED TO FIGURE OUT HOW TO CALCULATE FOR THE PHENO STUFF
// dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);

