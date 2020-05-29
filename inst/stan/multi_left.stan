functions {

#include /functions/cline_functions.stan

}

data {

#include /data/data_multinom.stan
#include /priors/priors_all.stan
#include /priors/priors_left_tail.stan
#include /priors/priors_f.stan

}

parameters{

#include /parameters/param_all.stan
#include /parameters/param_left_tail.stan
#include /parameters/param_f.stan

}

transformed parameters {
  vector<lower=0, upper=1>[N] p; // the allele frequency at each site
  vector[3] probs[N]; // The expected proportion of each genotype, based on p and f

  for ( i in 1:N )
  { // for each site
    if (transectDist[i] <= center - deltaL) { // if in left tail

      p[i] = left_tail(transectDist[i], center, width, pmin, pmax, deltaL, tauL, decrease);

    } else { // not in tail, use regular equation

      p[i] = gen_cline_eq(transectDist[i], center, width, pmin, pmax, decrease);

    }
  }

  for ( j in 1:N )
  { // take ps and fs and make genotype probs, which go into the likelihood
    probs[j,1] = p[j]^2 + f[j]*p[j]*(1-p[j]); // for AA
    probs[j,2] = 2*p[j]*(1-p[j])*(1-f[j]); // for Aa
    probs[j,3] = (1-p[j])^2 +f[j]*p[j]*(1-p[j]); // for aa
  }
}

model{
  // The statistical model

#include /model/model_f.stan
#include /model/model_left_tail.stan
#include /model/model_ps.stan
#include /model/model_width.stan
#include /model/model_center.stan

  if (ignoreData == 0) {
   for (i in 1:N ) { // for each site
  genos[i,] ~ multinomial(probs[i]); // the likelihood of genotype counts is multinomial for the 3 categories
  }
  }

}
generated quantities{

#include /generated_quantities/gen_quant_multinomial.stan

}
