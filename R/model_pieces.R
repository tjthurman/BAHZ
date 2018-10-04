
# Stan model pieces -------------------------------------------------------
# Already written Stan model code, for each of the models decribed. For each
# model, there's the code that comes before the priors, and the code that comes
# after it. The priors get added in the create_cline_models function, and output
# as model code to the global environment.


# Binomial, no tails, increasing ------------------------------------------
bi_none_inc_before_priors <- "// Stan model to fit a sigmoid cline of allele frequencies,
// with no introgression tails.
// The minimum and maximum alele frequencies at the tails of the cline
// are estimated from the data
data{
int<lower=1> N; // number of sites sampled
int nFocalAllele[N]; // number of alleles of the focal type, in this case hydara alleles
int nTotalAlleles[N]; // total number of alleles sampled (2*number of diploid individuals)
real transectDist[N]; // distance along transect
}
parameters{
real<lower=0> center; // the center of the cline, in km. Can't be negative
real<lower=0> width; // the width of the cline. Also can't be negative
real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
}
model{
// The statistical model
vector[N] p; // a vector of expected p values for each site\n"



bi_none_inc_after_priors <- "\nfor ( i in 1:N ) { // for each site
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

"


# Binomial, left tail, increasing -----------------------------------------
bi_left_inc_before_priors <- "// Stan model to fit a sigmoid cline of allele frequencies,
// with an introgression tail on the left
// The minimum and maximum alele frequencies at the tails of the cline
// are estimated from the data
data{
int<lower=1> N; // number of sites sampled
int nFocalAllele[N]; // number of alleles of the focal type, in this case hydara alleles
int nTotalAlleles[N]; // total number of alleles sampled (2*number of diploid individuals)
real transectDist[N]; // distance along transect
}
parameters{
real<lower=0> center; // the center of the cline, in km. Can't be negative
real<lower=0> width; // the width of the cline. Also can't be negative
real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
real<lower=0,upper=width> deltaL; // the location of the left tail.
real<lower=0,upper=1> tauL; // the ratio of the slope of the exponential tail to the slope of the sigmoid center. Must be between 0 and 1.
}
model{
// The statistical model
vector[N] p; // a vector of expected p values for each site\n"

bi_left_inc_after_priors <- "\nfor ( i in 1:N ) { // for each site
if(transectDist[i] <= center - deltaL) { // if the site is in the left tail
// predicted value of p is determined by the exponential equation for the left tail
p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width)));
}
else { // otherwise, the site is in the exponential center
// so use the regular sigmoid cline equation
p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
}
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
// with the expected p calculated above, and the same number of draws as the empirical data
for ( i in 1:N ) { // for each site
// calculate the predicted p for each site.
// These equations and if/else statments should be the same as above in the model
if (transectDist[i] <= center - deltaL) {
p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width)));
}
else {
p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
}
// and calculate the log likelihood of the observed allele counts, given the predicted p
log_lik[i] = binomial_lpmf(nFocalAllele[i] | nTotalAlleles[i], p[i]);
// and calculate an expected number of focal alleles, given the same number of samples
y_rep[i] = binomial_rng(nTotalAlleles[i], p[i]);
}
// calculate deviance of the model
dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);
}

"

# Binomial, right tail, increasing ----------------------------------------
bi_right_inc_before_priors <- "// Stan model to fit a sigmoid cline of allele frequencies,
// with an introgression tail on the right.
// The minimum and maximum alele frequencies at the tails of the cline
// are estimated from the data
data{
int<lower=1> N; // number of sites sampled
int nFocalAllele[N]; // number of alleles of the focal type, in this case hydara alleles
int nTotalAlleles[N]; // total number of alleles sampled (2*number of diploid individuals)
real transectDist[N]; // distance along transect
}
parameters{
real<lower=0> center; // the center of the cline, in km. Can't be negative
real<lower=0> width; // the width of the cline. Also can't be negative
real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
real<lower=0,upper=width> deltaR; // the location of the right tail.
real<lower=0,upper=1> tauR; // the ratio of the slope of the exponential tail to the slope of the sigmoid center. Must be between 0 and 1.
}
model{
// The statistical model
vector[N] p; // a vector of expected p values for each site\n"

bi_right_inc_after_priors <- "\nfor ( i in 1:N ) { // for each site
if(transectDist[i] >= center + deltaR) { // if the site is in the right tail
// predicted value of p is determined by the exponential equation for the right tail
p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
}
else { // otherwise, the site is in the exponential center
// so use the regular sigmoid cline equation
p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
}
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
// with the expected p calculated above, and the same number of draws as the empirical data
for ( i in 1:N ) { // for each site
// calculate the predicted p for each site.
// These equations and if/else statments should be the same as above in the model
if (transectDist[i] >= center + deltaR) {
p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
}
else {
p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
}
// and calculate the log likelihood of the observed allele counts, given the predicted p
log_lik[i] = binomial_lpmf(nFocalAllele[i] | nTotalAlleles[i], p[i]);
// and calculate an expected number of focal alleles, given the same number of samples
y_rep[i] = binomial_rng(nTotalAlleles[i], p[i]);
}
// calculate deviance of the model
dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);
}

"


# Binomial, mirror tails, increasing --------------------------------------
bi_mirror_inc_before_priors <- "// Stan model to fit a sigmoid cline of allele frequencies,
// with mirrored introgression tails.
// The minimum and maximum alele frequencies at the tails of the cline
// are estimated from the data
data{
int<lower=1> N; // number of sites sampled
int nFocalAllele[N]; // number of alleles of the focal type, in this case hydara alleles
int nTotalAlleles[N]; // total number of alleles sampled (2*number of diploid individuals)
real transectDist[N]; // distance along transect
}
parameters{
real<lower=0> center; // the center of the cline, in km. Can't be negative
real<lower=0> width; // the width of the cline. Also can't be negative
real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
real<lower=0,upper=width> deltaM; // the location of the tails.
real<lower=0,upper=1> tauM; // the ratio of the slope of the exponential tail to the slope of the sigmoid center. Must be between 0 and 1.
}
model{
// The statistical model
vector[N] p; // a vector of expected p values for each site\n"

bi_mirror_inc_after_priors <- "\nfor ( i in 1:N ) { // for each site
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
// with the expected p calculated above, and the same number of draws as the empirical data
for ( i in 1:N ) { // for each site
// calculate the predicted p for each site.
// These equations and if/else statments should be the same as above in the model
if (transectDist[i] <= center - deltaM) {
p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaM/width)))*exp((4*tauM*(transectDist[i] - center + deltaM)/width)/(1 + exp(-4*deltaM/width)));
}
else if (transectDist[i] >= center + deltaM) {
p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaM/width)))*exp((-4*tauM*(transectDist[i] - center - deltaM)/width)/(1 + exp(-4*deltaM/width))));
}
else {
p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
}
// and calculate the log likelihood of the observed allele counts, given the predicted p
log_lik[i] = binomial_lpmf(nFocalAllele[i] | nTotalAlleles[i], p[i]);
// and calculate an expected number of focal alleles, given the same number of samples
y_rep[i] = binomial_rng(nTotalAlleles[i], p[i]);
}
// calculate deviance of the model
dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);
}

"

# Binomial, independent tails, increasing ---------------------------------
bi_ind_inc_before_priors <- "// Stan model to fit a sigmoid cline of allele frequencies,
// with independent exponential tails
// The minimum and maximum alele frequencies at the tails of the cline
// are estimated from the data
data{
  int<lower=1> N; // number of sites sampled
  int nFocalAllele[N]; // number of alleles of the focal type, in this case hydara alleles
  int nTotalAlleles[N]; // total number of alleles sampled (2*number of diploid individuals)
  real transectDist[N]; // distance along transect
}
parameters{
  real<lower=0> center; // the center of the cline, in km. Can't be negative
  real<lower=0> width; // the width of the cline. Also can't be negative
  real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
  real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
  real<lower=0> deltaR; // the location of the right tail.
  real<lower=0,upper=1> tauR; // the ratio of the slope of the exponential tail to the slope of the sigmoid center. Must be between 0 and 1.
  real<lower=0,upper=width> deltaL; // same as deltaR, but for the left tail
  real<lower=0,upper=1> tauL; // the ratio of the slope of the exponential tail to the slope of the sigmoid center. Must be between 0 and 1.
}
model{
  // The statistical model
  vector[N] p; // a vector of expected p values for each site\n"

bi_ind_inc_after_priors <-"\nfor ( i in 1:N ) {
    // for each site
    if (transectDist[i] <= center - deltaL) { // if the site is in the left tail
      // predicted value of p is determined by the exponential equation for the left tail
      p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width)));
    }
    else if(transectDist[i] >= center + deltaR) { // if the site is in the right tail
      // predicted value of p is determined by the exponential equation for the right tail
      p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
    }
    else { // otherwise, the site is in the exponential center
      // so use the regular sigmoid cline equation
      p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
    }
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
  // with the expected p calculated above, and the same number of draws as the empirical data
  for ( i in 1:N ) {// for each site
    // calculate the predicted p for each site.
    // These equations and if/else statments should be the same as above in the model
    if (transectDist[i] <= center - deltaL) {
      p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width)));
    }
    else if (transectDist[i] >= center + deltaR) {
      p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
    }
    else {
      p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
    }
    // and calculate the log likelihood of the observed allele counts, given the predicted p
    log_lik[i] = binomial_lpmf(nFocalAllele[i] | nTotalAlleles[i], p[i]);
    // and calculate an expected number of focal alleles, given the same number of samples
    y_rep[i] = binomial_rng(nTotalAlleles[i], p[i]);
  }
  // calculate deviance of the model
  dev = (-2)*binomial_lpmf(nFocalAllele | nTotalAlleles, p);
}

"

# Multinomial, no tails, increasing ---------------------------------------
multi_none_inc_before_priors <- "// Stan model to fit a sigmoid cline of allele frequencies,
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
}
transformed parameters {
  vector<lower=0, upper=1>[N] p; // the allele frequency at each site
  vector[3] probs[N]; // The expected proportion of each genotype, based on p and f
  for (i in 1:N ) { // for each site
    p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
    probs[i,1] = p[i]^2 + f[i]*p[i]*(1-p[i]); // for AA
    probs[i,2] = 2*p[i]*(1-p[i])*(1-f[i]); // for Aa
    probs[i,3] = (1-p[i])^2 +f[i]*p[i]*(1-p[i]); // for aa
  }
}
model{
  // priors on parameters\n"

multi_none_inc_after_priors <- "\nfor (i in 1:N ) { // for each site
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

"


