// Stan model to fit a sigmoid cline of allele frequencies,
// with no introgression tails. 
// The minimum and maximum alele frequencies at the tails of the cline
// are estimated from the data
data{
  int<lower=1> N; // number of sites sampled
  int genos[N,3]; // number of individuals of each genotype, going AA, Aa, aa
  real transectDist[N]; // distance along transect
}
parameters{
  real<lower=0> center; // the center of the cline, in km. Can't be negative
  real<lower=0> width; // the width of the cline. Also can't be negative
  vector<lower=0, upper=1>[N] f; // inbreeding coefficieint for each site
  real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
  real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
  real<lower=0,upper=width> deltaM; // the location of the tails.
  real<lower=0,upper=1> tauM; // the ratio of the slope of the exponential tails to the slope of the sigmoid center. Must be between 0 and 1. 
}
transformed parameters {
  vector<lower=0, upper=1>[N] p; // the allele frequency at each site
  vector[3] probs[N]; // The expected proportion of each genotype, based on p and f
  for (i in 1:N ) { // for each site
    if (transectDist[i] <= center - deltaM) { // if the site is in the left tail
      // predicted value of p is determined by the exponential equation for the left tail
      p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaM/width)))*exp((4*tauM*(transectDist[i] - center + deltaM)/width)/(1 + exp(-4*deltaM/width)));
    }
    else if(transectDist[i] >= center + deltaM) { // if the site is in the right tail
      // predicted value of p is determined by the exponential equation for the right tail
      p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaM/width)))*exp((-4*tauM*(transectDist[i] - center - deltaM)/width)/(1 + exp(-4*deltaM/width))));
    } 
    else { // otherwise, the site is in the exponential center
      // so use the regular sigmoid cline equation
      p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
    }
    // Either way, the probs are calculated from the p and f values
      probs[i,1] = p[i]^2 + f[i]*p[i]*(1-p[i]); // for AA
      probs[i,2] = 2*p[i]*(1-p[i])*(1-f[i]); // for Aa
      probs[i,3] = (1-p[i])^2 +f[i]*p[i]*(1-p[i]); // for aa
  }
}
model{
  // priors on parameters
  f ~ uniform(0, 1); //prior for f. Uniform
  tauM ~ uniform(0,1); // prior for tauM
  deltaM ~ exponential(0.05); // prior for deltaM
  pmax ~ uniform(0.8 , 1); // prior for pmax
  pmin ~ uniform(0 , 0.2); // prior for pmin
  width ~ normal(50, 100); // prior for width. 
  center ~ normal(350, 150); //prior for center
  for (i in 1:N ) { // for each site
  genos[i,] ~ multinomial(probs[i]); // the likelihood of genotype counts is multinomial for the 3 categories
  }
}

generated quantities{
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
  // dev = (-2)*multinomial_lpmf(genos | probs);
}
