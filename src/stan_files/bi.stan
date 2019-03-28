#include /pre/license.stan

data {

#include /data/data_binom.stan
#include /priors/priors_all.stan
if (tails == 1) { //left
#include /priors/priors_left_tail.stan
} else if (tails == 2) { // right
#include /priors/priors_right_tail.stan
} else if (tails == 3) { //mirror
#include /priors/priors_mirror_tail.stan
}  else if (tails == 4) { //ind
#include /priors/priors_left_tail.stan
#include /priors/priors_right_tail.stan
} else { // no tails, do nothing.

}
}

parameters{

#include /parameters/param_all.stan
if (tails == 1) { //left
#include /parameters/param_left_tail.stan
} else if (tails == 2) { // right
#include /parameters/param_right_tail.stan
} else if (tails == 3) { //mirror
#include /parameters/param_mirror_tail.stan
}  else if (tails == 4) { //ind
#include /parameters/param_left_tail.stan
#include /parameters/param_right_tail.stan
} else { // no tails, do nothing.

}

}

transformed parameters {
  vector[N] p; // a vector of expected p values for each site
  for ( i in 1:N ) { // for each site
      if (decrease == 0) { // if increasing
        if (tails == 1) { // if left tail
          if (transectDist[i] <= center - deltaL) { // if in tail
             p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width)));
          } else { //regular equation.
              p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
          }
        } else if (tails == 2) { // if right tail
            if (transectDist[i] >= center + deltaR) { // if in tail
                p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
            } else { //regular equation.
                p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
          }
        } else if (tails == 3) { //if mirror tail
           if (transectDist[i] <= center - deltaM) { // if in left tail
             p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaM/width)))*exp((4*tauM*(transectDist[i] - center + deltaM)/width)/(1 + exp(-4*deltaM/width)));
           } else if (transectDist[i] >= center + deltaM) { // if in right tail
              p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaM/width)))*exp((-4*tauM*(transectDist[i] - center - deltaM)/width)/(1 + exp(-4*deltaM/width))));
        } else { // regular equation
          p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
        }
       } else if { // ind tail
        if (transectDist[i] <= center - deltaL) { // if in left tail
             p[i] =  pmin + (pmax-pmin)*(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width)));
           } else if (transectDist[i] >= center + deltaR) { // if in right tail
              p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
        } else { // regular equation
          p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
        }

      } else { // no tails, regular equation
        p[i] =   pmin + (pmax - pmin) * (exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width)));
      }

      }

    // // calculate the predicted p
    //  if (transectDist[i] <= center - deltaL) { // if the site is in the left tail
    //     // predicted value of p is determined by the exponential equation for the left tail
    //     }
    //     if (decrease == 1) {
    //        p[i] =  pmin + (pmax-pmin)*(1-(1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist[i] - center + deltaL)/width)/(1 + exp(-4*deltaL/width))));
    //     }
    //  } else if (transectDist[i] >= center + deltaR) { // if the site is in the right tail
    //     // predicted value of p is determined by the exponential equation for the right tail
    //     if (decrease == 1) {
    //       p[i] =  pmin + (pmax-pmin)*((1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist[i] - center - deltaR)/width)/(1 + exp(-4*deltaR/width))));
    //     }
    //  } else { // otherwise, the site is in the exponential center
    // // so use the regular sigmoid cline equation
    //       }
    //       if (decrease == 1) {
    //           p[i] =   pmin + (pmax - pmin) * (1-(exp(4*(transectDist[i] - center)/width)/(1 + exp(4 * (transectDist[i] - center)/width))));
    //       }
    //     }
    // }
}
}

model{
  // The statistical model

if (tails == 1) { //left
#include /model/model_left_tail.stan
} else if (tails == 2) { // right
#include /model/model_right_tail.stan
} else if (tails == 3) { //mirror
#include /model/model_mirror_tail.stan
}  else if (tails == 4) { //ind
#include /model/model_left_tail.stan
#include /model/model_right_tail.stan
} else { // no tails, do nothing.

}
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
