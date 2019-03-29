#include /pre/license.stan

data {

#include /data/data_multinom.stan
#include /priors/priors_all.stan
#include /priors/priors_left_tail.stan
#include /priors/priors_right_tail.stan
#include /priors/priors_mirror_tail.stan
#include /priors/priors_f.stan

}

parameters{

#include /parameters/param_all.stan
#include /parameters/param_left_tail.stan
#include /parameters/param_right_tail.stan
#include /parameters/param_mirror_tail.stan
#include /parameters/param_f.stan

}

transformed parameters {
  vector<lower=0, upper=1>[N] p; // the allele frequency at each site
  vector[3] probs[N]; // The expected proportion of each genotype, based on p and f
  for ( i in 1:N )
  { // for each site
      if (decrease == 0) // if increasing
      {
        if (tails == 1) // if left tail model
        {
          if (transectDist[i] <= center - deltaL) // if inside tail
          {
            p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width)));
          }
          else // not in tail, use regular equation.
          {
            p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
          }
        }
        else if (tails == 2) // if right tail model
        {
          if (transectDist[i] >= center + deltaR) // if inside tail
          {
            p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
          }
          else // not in tail, use regular equation.
          {
            p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
          }
        }
        else if (tails == 3) // if mirror tail model
        {
          if (transectDist[i] <= center - deltaM) // if inside left tail
          {
            p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaM/width)))*exp((4*tauM*(transectDist[i] - center + deltaM)/width)/(1 + exp(-4*deltaM/width)));
          }
          else if (transectDist[i] >= center + deltaM)  // if inside right tail
          {
            p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaM/width)))*exp((-4*tauM*(transectDist[i] - center - deltaM)/width)/(1 + exp(-4*deltaM/width))));
          }
          else
          {
            p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
          }
        }
        else if (tails == 4) // if ind tail model
        {
          if (transectDist[i] <= center - deltaL) // if in left tail
          {
            p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width)));
          }
          else if (transectDist[i] >= center + deltaR) // if in right tail
          {
            p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
          }
          else // not in tail, use regular equation.
          {
          p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
          }
        }
        else // no tail model, use regular equation
        {
          p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
       }
      }
       if (decrease == 1) // if decreasing # NEED TO MAKE ALL CHANGES
      {
        if (tails == 1) // if left tail model
        {
          if (transectDist[i] <= center - deltaL) // if inside tail
          {
            p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width))));
          }
          else // not in tail, use regular equation.
          {
            p[i] =   pmin + (pmax - pmin) * (1-(exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width))));
          }
        }
        else if (tails == 2) // if right tail model
        {
          if (transectDist[i] >= center + deltaR) // if inside tail
          {
            p[i] =  pmin + (pmax-pmin)*((1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
          }
          else // not in tail, use regular equation.
          {
            p[i] =   pmin + (pmax - pmin) * (1-(exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width))));
          }
        }
        else if (tails == 3) // if mirror tail model
        {
          if (transectDist[i] <= center - deltaM) // if inside left tail
          {
            p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaM/width)))*exp((4*tauM*(transectDist[i] - center + deltaM)/width)/(1 + exp(-4*deltaM/width))));
          }
          else if (transectDist[i] >= center + deltaM)  // if inside right tail
          {
            p[i] =  pmin + (pmax-pmin)*((1/(1 + exp(4*deltaM/width)))*exp((-4*tauM*(transectDist[i] - center - deltaM)/width)/(1 + exp(-4*deltaM/width))));
          }
          else
          {
            p[i] =   pmin + (pmax - pmin) * (1-(exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width))));
          }
        }
        else if (tails == 4) // if ind tail model
        {
          if (transectDist[i] <= center - deltaL) // if in left tail
          {
            p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width))));
          }
          else if (transectDist[i] >= center + deltaR) // if in right tail
          {
            p[i] =  pmin + (pmax-pmin)*((1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
          }
          else // not in tail, use regular equation.
          {
          p[i] =   pmin + (pmax - pmin) * (1-(exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width))));
          }
        }
        else // no tail model, use regular equation
        {
          p[i] =   pmin + (pmax - pmin) * (1-(exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width))));
       }
      }
  } // That all calculates p
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
#include /model/model_mirror_tail.stan
#include /model/model_left_tail.stan
#include /model/model_right_tail.stan
#include /model/model_ps.stan
#include /model/model_width.stan
#include /model/model_center.stan

  for (i in 1:N ) { // for each site
  genos[i,] ~ multinomial(probs[i]); // the likelihood of genotype counts is multinomial for the 3 categories
  }
}
generated quantities{

#include /generated_quantities/gen_quant_multinomial.stan

}
