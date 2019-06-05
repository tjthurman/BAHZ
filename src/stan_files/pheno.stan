#include /pre/license.stan

data {

#include /data/data_pheno.stan
#include /priors/priors_all.stan

}

parameters{

#include /parameters/param_pheno.stan

}

transformed parameters{
vector[K] p; // the expected phenotype for each site
for (i in 1:K) {
  p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
  }
}


model{
  // The statistical model

sigma ~ normal(0, 100);
#include /model/model_ps.stan
#include /model/model_width.stan
#include /model/model_center.stan

for (k in 1:K) {
  pheno[starts[k]:ends[k]] ~ normal(p[k], sigma);
  }
}

generated quantities{

#include /generated_quantities/gen_quant_pheno.stan

}
