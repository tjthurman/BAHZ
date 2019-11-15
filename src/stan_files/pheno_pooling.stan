#include /pre/license.stan

data {

#include /data/data_pheno.stan
#include /priors/priors_all.stan

}

parameters{

#include /parameters/param_pheno.stan
#include /parameters/param_pheno_pooling.stan

}

transformed parameters{
// expected phenotype per site
vector[K] p; // the expected phenotype for each site
for (i in 1:K) {
  if (decrease == 0) // if increasing
      {
        p[i] = pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
        }
  if (decrease == 1) // if decreasing
      {
        p[i] =   pmin + (pmax - pmin) * (1-(exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width))));
        }
}
}

model{
  // The statistical model

sigma_sigma ~ cauchy(0, 50);
global_sigma ~ cauchy(0, 50);
site_sigma ~ normal(global_sigma, sigma_sigma);
#include /model/model_ps.stan
#include /model/model_width.stan
#include /model/model_center.stan

for (k in 1:K) {
  pheno[starts[k]:ends[k]] ~ normal(p[k], site_sigma[k]);
  }
}

generated quantities{

#include /generated_quantities/gen_quant_pheno_not_constant.stan

}
