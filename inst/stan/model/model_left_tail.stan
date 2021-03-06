if (p_dist_tauL == 1) { // uniform
  tauL ~ uniform(p_tauL_1, p_tauL_2);
}
if (p_dist_tauL == 3) { // beta
  tauL ~ beta(p_tauL_1, p_tauL_2);
}
if (p_dist_deltaL == 1) { // uniform
  deltaL ~ uniform(p_deltaL_1, p_deltaL_2);
}
if (p_dist_deltaL == 2) { // exponential
  deltaL ~ exponential(p_deltaL_1);
}
