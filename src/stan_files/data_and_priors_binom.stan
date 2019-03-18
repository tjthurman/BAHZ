data{
  int<lower=1> N; // number of sites sampled
  int nFocalAllele[N]; // number of alleles of the focal type, in this case hydara alleles
  int nTotalAlleles[N]; // total number of alleles sampled (2*number of diploid individuals)
  real transectDist[N]; // distance along transect
  real p_m_center; // the mean of the prior distribution for center
  real p_sd_center; // the sd of the prior distribution for center
  real p_m_width; // the mean of the prior distribution for width
  real p_sd_width; // the sd of the prior distribution for width
  real p_scale_width; // the sd of the prior distribution for width
  real p_l_min; // the lower bound of the prior distribution for pmin
  real p_u_min; // the upper bound of the prior distribution for pmin
  real p_l_max; // the lower bound of the prior distribution for pmin
  real p_u_max; // the upper bound of the prior distribution for pmin
}
