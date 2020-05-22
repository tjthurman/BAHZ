real gen_cline_eq(real dist, real c, real w, real pmi, real pmx, int dec) {
  real int_res;
  real final_res;
  int_res = exp(4*(dist - c)/w)/(1 + exp(4 * (dist - c)/w));
  if (dec == 0) {
    final_res = pmi + (pmx - pmi) * int_res;
  }
  if (dec == 1) {
    final_res = pmi + (pmx - pmi) * (1-int_res);
  }
  return final_res;
}

real left_tail(real dist, real c, real w, real pmi, real pmx, real dl, real tl, int dec) {
  real int_res;
  real final_res;

  int_res = (1/(1 + exp(4*dl/w)))*exp((4*tl*(dist - c + dl)/w)/(1 + exp(-4*dl/w)));

  if (dec == 0) {
    final_res = pmi + (pmx - pmi) * int_res;
  }
  if (dec == 1) {
    final_res = pmi + (pmx - pmi) * (1-int_res);
  }
  return final_res;
}

real right_tail(real dist, real c, real w, real pmi, real pmx, real dr, real tr, int dec) {
  real int_res;
  real final_res;

  int_res = (1-(1/(1 + exp(4*dr/w)))*exp((-4*tr*(dist - c - dr)/w)/(1 + exp(-4*dr/w))));

  if (dec == 0) {
    final_res = pmi + (pmx - pmi) * int_res;
  }
  if (dec == 1) {
    final_res = pmi + (pmx - pmi) * (1-int_res);
  }
  return final_res;
}


