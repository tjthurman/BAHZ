int<lower=1> N; // number of sites sampled
int nFocalAllele[N]; // number of alleles of the focal type, in this case hydara alleles
int nTotalAlleles[N]; // total number of alleles sampled (2*number of diploid individuals)
real transectDist[N]; // distance along transect
int<lower = 0, upper = 1> decrease; // is cline ascending or descending (guessed from data when passed to Stan)
int<lower = 0, upper = 4> tails; // Which tail model to use? 0 = none, 1 = left, 2 = right, 3 = mirror, 4 = ind
