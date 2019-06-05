if (p_dist_pmin == 0) { // normal
  pmin ~ normal(p_min_1, p_min_2);
}
if (p_dist_pmin == 1) { // uniform
  pmin ~ uniform(p_min_1, p_min_2);
}
if (p_dist_pmin == 3) { // beta
  pmin ~ beta(p_min_1, p_min_2);
}
if (p_dist_pmax == 0) { // normal
  pmax ~ normal(p_max_1, p_max_2);
}
if (p_dist_pmax == 1) { // uniform
  pmax ~ uniform(p_max_1, p_max_2);
}
if (p_dist_pmax == 3) { // beta
  pmax ~ beta(p_max_1, p_max_2);
}
