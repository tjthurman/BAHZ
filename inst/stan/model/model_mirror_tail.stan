if (p_dist_tauM == 1) { // uniform
  tauM ~ uniform(p_tauM_1, p_tauM_2);
}
if (p_dist_tauM == 3) { // beta
  tauM ~ beta(p_tauM_1, p_tauM_2);
}
deltaM ~ exponential(p_deltaM_1); // prior for deltaM
