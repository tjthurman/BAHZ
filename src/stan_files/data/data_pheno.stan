int<lower=0> N; // number of individuals
int<lower=0> K; // number of sites
vector[N] pheno; // phenotype for an individual
int s[K]; //number of individuals sampled per site
real transectDist[K]; // distance along transect for a site
int starts[K]; // the indexes of the total # of individuals where each site starts
int ends[K]; // the indexes of the total # of individuals where each site ends
