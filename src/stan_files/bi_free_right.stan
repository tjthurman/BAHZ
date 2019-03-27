#include /pre/license.stan

data {

#include /data/data_binom.stan
#include /priors/priors_all.stan
#include /priors/priors_right_tail.stan

}

parameters{

#include /parameters/param_all.stan
#include /parameters/param_right_tail.stan

}

transformed parameters {
  vector[N] p; // a vector of expected p values for each site
  for ( i in 1:N ) { // for each site
    // calculate the predicted p
     if (transectDist[i] >= center + deltaR) { // if the site is in the right tail
        // predicted value of p is determined by the exponential equation for the right tail
        if (decrease == 0) {
          p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
        }
        if (decrease == 1) {
          p[i] =  pmin + (pmax-pmin)*((1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
        }
     } else { // otherwise, the site is in the exponential center
    // so use the regular sigmoid cline equation
          if (decrease == 0) {
              p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
          }
          if (decrease == 1) {
              p[i] =   pmin + (pmax - pmin) * (1-(exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width))));
          }
        }
    }
}

model{
  // The statistical model

#include /model/model_right_tail.stan
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
