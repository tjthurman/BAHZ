if (p_dist_tauR == 1) { // uniform
  tauR ~ uniform(p_tauR_1, p_tauR_2);
}
if (p_dist_tauR == 3) { // beta
  tauR ~ beta(p_tauR_1, p_tauR_2);
}
if (p_dist_deltaR == 1) { // uniform
  deltaR ~ uniform(p_deltaR_1, p_deltaR_2);
}
if (p_dist_deltaR == 2) { // exponential
  deltaR ~ exponential(p_deltaR_1);
}
