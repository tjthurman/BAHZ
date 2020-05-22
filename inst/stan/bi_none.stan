functions {

#include /functions/cline_functions.stan

}

data {

#include /data/data_binom.stan
#include /priors/priors_all.stan

}

parameters{

#include /parameters/param_all.stan

}

transformed parameters {
  vector<lower=0, upper=1>[N] p; // the allele frequency at each site
  for ( i in 1:N )
  { // for each site
    p[i] = gen_cline_eq(transectDist[i], center, width, pmin, pmax, decrease);
  }
}

model{
  // The statistical model


#include /model/model_ps.stan
#include /model/model_width.stan
#include /model/model_center.stan


  // and the likelihood: observed allele counts follow a binomial liklihood,
  // based on the number of alleles sampled and the estimated allele frequency.
  nFocalAllele ~ binomial(nTotalAlleles, p);
}
generated quantities{

#include /generated_quantities/gen_quant_binomial.stan

}
