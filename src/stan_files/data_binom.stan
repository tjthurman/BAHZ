data{
  int<lower=1> N; // number of sites sampled
  int nFocalAllele[N]; // number of alleles of the focal type, in this case hydara alleles
  int nTotalAlleles[N]; // total number of alleles sampled (2*number of diploid individuals)
  real transectDist[N]; // distance along transect
}
