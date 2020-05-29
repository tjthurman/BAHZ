if (p_dist_center == 0) { // normal
  center ~ normal(p_center_1, p_center_2); //prior for center
}
if (p_dist_center == 1) { // uniform
  center ~ uniform(p_center_1, p_center_2); //prior for center
}
