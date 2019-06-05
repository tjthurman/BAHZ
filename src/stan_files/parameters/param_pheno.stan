real center; // the center of the cline, in km.
real<lower=0> width; // the width of the cline. Can't be negative
real pmin; // minimum mean pheno value
real pmax; // maximum mean pheno value
real<lower =0> sigma; // phenotypic variance (constant for now, must be positive)
