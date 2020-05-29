int<lower=0> N; // number of individuals
int<lower=0> K; // number of sites
real pheno[N]; // phenotype for an individual
real indDist[N]; // distance for an individual
int<lower = 0, upper = 1> decrease; // is cline ascending or descending (guessed from data when passed to Stan)
int s[K]; //number of individuals sampled per site
real transectDist[K]; // distance along transect for a site
int starts[K]; // the indexes of the total # of individuals where each site starts
int ends[K]; // the indexes of the total # of individuals where each site ends
int<lower = 0, upper = 1> ignoreData; // should we ignore the data to do prior predictive modelling?
