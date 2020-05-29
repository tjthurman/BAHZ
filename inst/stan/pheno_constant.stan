functions {

#include /functions/cline_functions.stan

}

data {

#include /data/data_pheno.stan
#include /priors/priors_all.stan

}

parameters{

#include /parameters/param_pheno.stan
#include /parameters/param_pheno_constant.stan

}

transformed parameters{
  // expected phenotype per site
  vector[K] p; // the expected phenotype for each site
  for ( i in 1:K )
  { // for each site
    p[i] = gen_cline_eq(transectDist[i], center, width, pmin, pmax, decrease);
  }
}

model{
  // The statistical model

global_sigma ~ cauchy(0, 50);
#include /model/model_ps.stan
#include /model/model_width.stan
#include /model/model_center.stan

if (ignoreData == 0) {
for (k in 1:K) {
  pheno[starts[k]:ends[k]] ~ normal(p[k], global_sigma);
  }
}

}

generated quantities{

#include /generated_quantities/gen_quant_pheno_constant.stan

}
