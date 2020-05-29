if (p_dist_f == 1) { //uniform
  f ~ uniform(p_f_1, p_f_2); // prior for f
}
if (p_dist_f == 3) { //beta
  f ~ beta(p_f_1, p_f_2); // prior for f
}
