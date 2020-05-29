int<lower=1> N; // number of sites sampled
int genos[N,3]; // number of individuals of each genotype, going AA, Aa, aa
real transectDist[N]; // distance along transect
int<lower = 0, upper = 1> decrease; // is cline ascending or descending (guessed from data when passed to Stan)
int<lower = 0, upper = 1> ignoreData; // should we ignore the data to do prior predictive modelling?


