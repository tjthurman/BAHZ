context("Simulate geno data")

# Test error handling
test_that("sim_geno_cline checks ind and fis for type (vectors)", {
  expect_error(sim_geno_cline(transect_distances = seq(0,600,20),
                              n_ind = 40,
                              Fis = as.matrix(40.5, 30),
                              decrease = F, center = 300, width = 80), "vector")
  expect_error(sim_geno_cline(transect_distances = seq(0,600,20),
                              n_ind = as.matrix(40.5, 30),
                              Fis = 0,
                              decrease = F, center = 300, width = 80), "vector")

})
test_that("sim_geno_cline checks ind and fis for being numeric", {
  expect_error(sim_geno_cline(transect_distances = seq(0,600,20),
                              n_ind = 40,
                              Fis = "A",
                              decrease = F, center = 300, width = 80), "Fis must be numeric")
  expect_error(sim_geno_cline(transect_distances = seq(0,600,20),
                              n_ind = "A",
                              Fis = 0,
                              decrease = F, center = 300, width = 80), "n_ind must be numeric")

})
test_that("sim_geno_cline checks ind and fis for proper length", {
  expect_error(sim_geno_cline(transect_distances = seq(0,600,20),
                              n_ind = 40,
                              Fis = c(40, 30),
                              decrease = F, center = 300, width = 80), "length")
  expect_error(sim_geno_cline(transect_distances = seq(0,600,20),
                              n_ind = c(40, 30),
                              Fis = 0,
                              decrease = F, center = 300, width = 80), "length")

})
test_that("sim_geno_cline checks that ind and fis are in proper range", {
  expect_error(sim_geno_cline(transect_distances = seq(0,600,20),
                              n_ind = 20,
                              Fis = -0.1,
                              decrease = F, center = 300, width = 80), "less")
  expect_error(sim_geno_cline(transect_distances = seq(0,600,20),
                              n_ind = 20,
                              Fis = 1.5,
                              decrease = F, center = 300, width = 80), "greater")
  expect_error(sim_geno_cline(transect_distances = seq(0,600,20),
                              n_ind = 0,
                              Fis = 0.5,
                              decrease = F, center = 300, width = 80), "less")
})

# Test cline generation

# First, generate a bunch of clines to be used a references later,
# which test the full functionality of cline generation.
# Will save those as a reference object to load, and comment out the
# saving step.
# For now, will do a bunch of versions of the most complicated
# clines possible.
set.seed(10)
a <-sim_geno_cline(transect_distances = seq(0,600,20), n_ind = 40, Fis = 0,
                   decrease = F,
                   center = 300, width = 80,
                   pmin = 0.05,
                   pmax = 0.95,
                   deltaL = 10, tauL = 0.5, deltaR = 11, tauR = 0.6)
b <-sim_geno_cline(transect_distances = seq(0,600,20), n_ind = 40, Fis = 0,
                   decrease = T,
                   center = 300, width = 80,
                   pmin = 0.05,
                   pmax = 0.95,
                   deltaL = 10, tauL = 0.5, deltaR = 11, tauR = 0.6)
c <-sim_geno_cline(transect_distances = seq(-300,300,20), n_ind = 40, Fis = 0,
                   decrease = F,
                   center = 300, width = 80,
                   pmin = 0.05,
                   pmax = 0.95,
                   deltaL = 10, tauL = 0.5, deltaR = 11, tauR = 0.6)
d <-sim_geno_cline(transect_distances = seq(-300,300,20), n_ind = 40, Fis = 0,
                   decrease = T,
                   center = 300, width = 80,
                   pmin = 0.05,
                   pmax = 0.95,
                   deltaL = 10, tauL = 0.5, deltaR = 11, tauR = 0.6)
e <-sim_geno_cline(transect_distances = seq(-300,300,20),
                   n_ind = rnorm(length(seq(-300,300,20)), 50, 2),
                   Fis = runif(length(seq(-300,300,20)), 0, 1),
                   decrease = F,
                   center = 300, width = 80,
                   pmin = 0.05,
                   pmax = 0.95,
                   deltaL = 10, tauL = 0.5, deltaR = 11, tauR = 0.6)
f <-sim_geno_cline(transect_distances = seq(-300,300,20),
                   n_ind = rnorm(length(seq(-300,300,20)), 50, 2),
                   Fis = runif(length(seq(-300,300,20)), 0, 1),
                   decrease = T,
                   center = 300, width = 80,
                   pmin = 0.05,
                   pmax = 0.95,
                   deltaL = 10, tauL = 0.5, deltaR = 11, tauR = 0.6)
# for (object in c(expr(a), expr(b),expr(c),expr(d),expr(e),expr(f))) {
#  assign(x = paste(object, ".ref", sep = ""), value = eval(object))
# }
# save(a.ref,b.ref,c.ref,d.ref,e.ref,f.ref, file = "tests/testthat/ref_simulated_geno_clines.Rda")
load("ref_simulated_geno_clines.Rda")

test_that("clines generated by sim_geno_cline are correct", {
  expect_equal(a, a.ref)
  expect_equal(b, b.ref)
  expect_equal(c, c.ref)
  expect_equal(d, d.ref)
  expect_equal(e, e.ref)
  expect_equal(f, f.ref)
})
