if (p_dist_width == 0) { // normal
  width ~ normal(p_width_1, p_width_2); //prior for width
}
if (p_dist_width == 1) { // uniform
  width ~ uniform(p_width_1, p_width_2); //prior for width
}
