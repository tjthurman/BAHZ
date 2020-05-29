if (p_dist_tauM == 1) { // uniform
  tauM ~ uniform(p_tauM_1, p_tauM_2);
}
if (p_dist_tauM == 3) { // beta
  tauM ~ beta(p_tauM_1, p_tauM_2);
}
if (p_dist_deltaM == 1) { // uniform
  deltaM ~ uniform(p_deltaM_1, p_deltaM_2);
}
if (p_dist_deltaM == 2) { // exponential
  deltaM ~ exponential(p_deltaM_1);
}
