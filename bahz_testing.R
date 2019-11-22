
# Load packages -----------------------------------------------------------
rm(list = ls())
library(devtools)
detach("package:bahz", unload=TRUE)
install_github("tjthurman/BAHZ")
install.packages("~/Documents/Work/PhD:McGill/Projects/BAHZ/", repos = NULL, type = "source")
library(bahz)
library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(loo)
# Getting the minimal pre-compiled model to run -------------------
# Generate a dataset
data <- sim_geno_cline(transect_distances = seq(0,500,20), n_ind = 30, Fis = 0,
                    decrease = T, center = 238, width = 66, pmin = 0.03, pmax = .95, tauL = 0.5, deltaL = 12)

set.seed(22)
data <- sim_geno_cline(transect_distances = seq(-300,300,20), n_ind = 40, Fis = 0,
                       decrease = F, center = 10, width = 35, pmin = 0.08, pmax = .95)

plot(x = data$transectDist, y = data$emp.p)
lines(x = data$transectDist, y = data$cline.p)
data2 <- rbind(data[1,])

make_prior_config()
# Fit the model, binomial
fit_none_b <- fit_geno_cline(data = data, prior_file = "all_betas.yaml",
                           type = "bi", tails = "none")
fit_left_b <- fit_geno_cline(data = data, prior_file = "all_betas.yaml",
                      type = "bi", tails = "left")
fit_right_b <- fit_geno_cline(data = data, prior_file = "all_betas.yaml",
                           type = "bi", tails = "right")
fit_mirror_b <- fit_geno_cline(data = data, prior_file = "all_betas.yaml",
                           type = "bi", tails = "mirror")
fit_ind_b <- fit_geno_cline(data = data, prior_file = "all_betas.yaml",
                            type = "bi", tails = "ind")

plot_geno_cline(fit_none_b, data = data, add.obs.freqs = T, col = "red")
plot_geno_cline(fit_none_b, data = data, add.obs.freqs = T, col = "xxx")
plot_geno_cline(fit_none_b, data = dplyr::slice(data, -1), add.obs.freqs = T, point.col = "red")


plot_geno_cline(fit_left_b, data = data, add.obs.freqs = T, col = "red")
plot_geno_cline(fit_right_b, data = data, add.obs.freqs = T, col = "red")
plot_geno_cline(fit_mirror_b, data = data, add.obs.freqs = T, col = "red")
plot_geno_cline(fit_ind_b, data = data, add.obs.freqs = T, col = "red")

plot_geno_cline(fit_none_b, data = data, main = "test", col = "red", xlab = "distance", ylab = "allele frequency")

general_cline_eqn(c(8,10), decrease = F, center = 5, width = 8)

cline_summary(fit_none_b)
cline_summary(fit_left_b)
cline_summary(fit_right_b)
cline_summary(fit_mirror_b)
cline_summary(fit_ind_b)

z <- fit_geno_cline(data = data, prior_file = "prior_config_test1.yaml", type = "bi", tails = "none", chains = 1)



inits <- prep_init_list("prior_config_test1.yaml", tails = "none", chains = as.integer(1), type = "geno")
priors <- prep_prior_list("prior_config_test1.yaml")

z@inits


plot(predict_geno_cline(fit_none_b, data = data)$transectDist, predict_geno_cline(fit_none_b, data = data)$p, type = "l")


lines(-300:300, plot_cline(fit_left_b)$p, type = "l", col = "red")
lines(-300:300, plot_cline(fit_right_b)$p, type = "l", col = "blue")
lines(-300:300, plot_cline(fit_mirror_b)$p, type = "l", col = "orange")
lines(-300:300, plot_cline(fit_ind_b)$p, type = "l", col = "green")

# Fit the model, binomial
fit_none_m <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                           type = "bi", tails = "none")
fit_left_m <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                           type = "multi", tails = "left")
fit_right_m <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                            type = "multi", tails = "right")
fit_mirror_m <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                             type = "multi", tails = "mirror")
fit_ind_m <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                          type = "multi", tails = "ind")

z <- predict_cline(stanfit = fit_none_m, distance = 0:500, confidence = T, prob = 0.94)
z2 <- predict_cline(stanfit = fit_left_m, distance = 0:500, confidence = T, prob = 0.95)
z3 <- predict_cline(stanfit = fit_right_m, distance = 0:500, confidence = T, prob = 0.95)
z4 <- predict_cline(stanfit = fit_mirror_m, distance = 0:500, confidence = T, prob = 0.95)
z5 <- predict_cline(stanfit = fit_ind_m, distance = 0:500, confidence = T, prob = 0.95)


install.packages("memoise")

library(memoise)

mem_pred_cline <- memoise(predict_cline)

test <- memoise(function(x) {x + 1})

test(1)

z <- predict_cline(stanfit = fit_none_m, distance = -300:300, confidence = T, prob = 0.95)
z2 <- predict_cline(stanfit = fit_none_m, distance = -300:300, confidence = T, prob = 0.95)

data$transectDist[2] <- -281

plot_geno_cline(stanfit = fit_none_m, data = data, add.obs.freqs = T, confidence = T)
plot_geno_cline(stanfit = fit_none_m, data = data, add.obs.freqs = T, confidence = T,
                cline.col = "orange", point.col = "yellow", prob = 0.5,
                lwd = 4, pch = 22, cex = 4, bg = "green")

par()

ggplot() +
  geom_ribbon(fill = "grey90",
              aes(x = transectDist, ymin = low_0.95_HPDI, ymax = up_0.95_HPDI),
              data = z) +
  geom_point(aes(x = transectDist, y = emp.p), data = data) +
  geom_line(aes(x = transectDist, y = p), data = z) +
  ggtitle("none")

ggplot() +
  geom_ribbon(fill = "grey90",
              aes(x = transectDist, ymin = low_0.95_HPDI, ymax = up_0.95_HPDI),
              data = z2) +
  geom_point(aes(x = transectDist, y = emp.p), data = data) +
  geom_line(aes(x = transectDist, y = p), data = z2) +
  ggtitle("left")

ggplot() +
  geom_ribbon(fill = "grey90",
              aes(x = transectDist, ymin = low_0.95_HPDI, ymax = up_0.95_HPDI),
              data = z3) +
  geom_point(aes(x = transectDist, y = emp.p), data = data) +
  geom_line(aes(x = transectDist, y = p), data = z3) +
  ggtitle("right")

ggplot() +
  geom_ribbon(fill = "grey90",
              aes(x = transectDist, ymin = low_0.95_HPDI, ymax = up_0.95_HPDI),
              data = z4) +
  geom_point(aes(x = transectDist, y = emp.p), data = data) +
  geom_line(aes(x = transectDist, y = p), data = z4) +
  ggtitle("mirror")

ggplot() +
  geom_ribbon(fill = "grey90",
              aes(x = transectDist, ymin = low_0.95_HPDI, ymax = up_0.95_HPDI),
              data = z5) +
  geom_point(aes(x = transectDist, y = emp.p), data = data) +
  geom_line(aes(x = transectDist, y = p), data = z5) +
  ggtitle("ind")

rethinking::compare(fit)

cline_summary(fit_left_m)

prep_init_list("prior_config_template.yaml", tails = "ind", chains = as.integer(1))

fit_none_m@inits
fit_left_m@inits

cline_summary(fit_none_m)
cline_summary(fit_left_m)
cline_summary(fit_right_m)
cline_summary(fit_mirror_m)
cline_summary(fit_ind_m)


xs <- seq(-300,300,20)
ys <- rep(NA, length(xs))
for (i in 1:length(xs)) {
  ys[i] <- general_cline_eqn(transectDist = xs[i], decrease = T,
                         center = 7.49, width = 40.20, pmin = 0.03, pmax = 0.96,
                        deltaL = 0.2, tauL = 0.65, deltaR = 0.2, tauR = 0.65)
}
plot(x = data$transectDist, y = data$emp.p)
lines(x = data$transectDist, y = data$cline.p)
lines(x = xs, y = ys, col = "red")

plot(z1)

z1 <- loo::loo(fit_none_m, r_eff = relative_eff(fit_none_m))
z2 <- loo::loo(fit_left_m, r_eff = relative_eff(fit_left_m))
z3 <- loo::loo(fit_right_m, r_eff = relative_eff(fit_right_m))
z4 <- loo::loo(fit_mirror_m, r_eff = relative_eff(fit_mirror_m))
z5 <- loo::loo(fit_ind_m, r_eff = relative_eff(fit_ind_m))

rethinking::compare(fit_none, fit_left, fit_right, fit_mirror, fit_ind)
loo::compare(z1, z2, z3, z4, z5)

# testing the new width check in prep_init_list
i <- 8000
z <- prep_init_list("prior_config_template.yaml", tails = "right", chains = as.integer(i))
widths <- rep(as.numeric(NA), times = i)
for (width in 1:i) {
  widths[width] <- z[[width]]$width
}
unique(unlist(widths) < 0) # Never makes a value below 0.

test <- function(...) {
  stan.args <- names(sapply(match.call(), deparse))[-1]
  if (c("control") %in% stan.args == F) {
    result <- list(adapt_delta = 0.95)
  } else{
    result <- which(stan.args == "control")
  }
}


# Will stan run with no priors? -------------------------------------------

# Yes, stan will run with no priors! I think the default, then, is that it
# puts uniform -inf, Inf priors, which against the whole point of using stan,
# really. Still, good to test.


# Re-scaling of width? ----------------------------------------------------

# Is it possible to rescale the width parameter as well, using 4/width

# Testing on how to convert from normal distribution on width scale
# to normal distribution on w scale.

# After a lot of testing, the answer seems to be no.
# It would be nice to re-scale width to ease the starting and allow it
# to start from default values. But, I'm nt sure there's any way to convert
# the user-supplied prior to a sensible prior on the 4/width scale.


