real center; // the center of the cline, in km.
real<lower=0> width; // the width of the cline. Also can't be negative
real<lower=0,upper=1> pmin; // minimum allele frequency. must be between 0 and 1
real<lower=0,upper=1> pmax; // maximum allele frequency. must be between 0 and 1
