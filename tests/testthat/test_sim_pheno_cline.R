context("Simulate pheno data")

# Test error handling
test_that("sim_pheno_cline checks ind and sigma for type (vectors)", {
  expect_error(sim_pheno_cline(transect_distances = seq(0,600,20),
                              n_ind = 40,
                              sigma = as.matrix(40.5, 30),
                              decrease = F, center = 300, width = 80, pmin = 10,
                              pmax = 20), "vector")
  expect_error(sim_pheno_cline(transect_distances = seq(0,600,20),
                              n_ind = as.matrix(40.5, 30),
                              sigma = 0,
                              decrease = F, center = 300, width = 80,
                              pmin = 10, pmax = 20), "vector")

})
test_that("sim_pheno_cline checks ind and sigma for being numeric", {
  expect_error(sim_pheno_cline(transect_distances = seq(0,600,20),
                              n_ind = 40,
                              sigma = "A",
                              decrease = F, center = 300, width = 80,
                              pmin = 10, pmax = 20), "sigma must be numeric")
  expect_error(sim_pheno_cline(transect_distances = seq(0,600,20),
                              n_ind = "A",
                              sigma = 0,
                              decrease = F, center = 300, width = 80,
                              pmin = 10, pmax = 20), "n_ind must be numeric")

})
test_that("sim_pheno_cline checks ind and sigma for proper length", {
  expect_error(sim_pheno_cline(transect_distances = seq(0,600,20),
                              n_ind = 40,
                              sigma = c(40, 30),
                              decrease = F, center = 300, width = 80,
                              pmin = 10, pmax = 20), "length")
  expect_error(sim_pheno_cline(transect_distances = seq(0,600,20),
                              n_ind = c(40, 30),
                              sigma = 0,
                              decrease = F, center = 300, width = 80,
                              pmin = 10, pmax = 20), "length")

})
test_that("sim_pheno_cline checks that ind and sigma are in proper range", {
  expect_error(sim_pheno_cline(transect_distances = seq(0,600,20),
                              n_ind = 20,
                              sigma = -0.1,
                              decrease = F, center = 300, width = 80,
                              pmin = 10, pmax = 20), "less")
  expect_error(sim_pheno_cline(transect_distances = seq(0,600,20),
                              n_ind = 0,
                              sigma = 0.5,
                              decrease = F, center = 300, width = 80,
                              pmin = 10, pmax = 20), "less")
  expect_error(sim_pheno_cline(transect_distances = seq(0,600,20),
                               n_ind = 0,
                               sigma = 0.5,
                               decrease = F, center = 300, width = 80,
                               pmin = 20, pmax = 10), "less")
})


# Test cline generation

# First, generate a bunch of clines to be used a references later,
# which test the full functionality of cline generation.
# Will save those as a reference object to load, and comment out the
# saving step.
# For now, will do a bunch of versions of the most complicated
# clines possible.
set.seed(10)
a <-sim_pheno_cline(transect_distances = seq(0,600,20), n_ind = 40, sigma = 8,
                   decrease = F,
                   center = 300, width = 80,
                   pmin = 10, pmax = 20)
b <-sim_pheno_cline(transect_distances = seq(0,600,20), n_ind = 40, sigma = 12,
                    decrease = F,
                    center = 300, width = 10,
                    pmin = 10, pmax = 20)
c <-sim_pheno_cline(transect_distances = seq(-300,300,20), n_ind = 40, sigma = 8,
                   decrease = T,
                   center = 300, width = 80,
                   pmin = -10,
                   pmax = 10)
d <-sim_pheno_cline(transect_distances = seq(-300,300,20), n_ind = 40, sigma = 30,
                   decrease = T,
                   center = 300, width = 80,
                   pmin = 10,
                   pmax = 12)
e <-sim_pheno_cline(transect_distances = seq(-300,300,20),
                   n_ind = as.integer(rnorm(length(seq(-300,300,20)), 50, 2)),
                   sigma = runif(length(seq(-300,300,20)), 4, 8),
                   decrease = F,
                   center = 300, width = 80,
                   pmin = 20,
                   pmax = 30)
f <-sim_pheno_cline(transect_distances = seq(-300,300,20),
                   n_ind = as.integer(rnorm(length(seq(-300,300,20)), 50, 2)),
                   sigma = runif(length(seq(-300,300,20)), 30, 31),
                   decrease = T,
                   center = 300, width = 80,
                   pmin = 100,
                   pmax = 200)
# for (object in c(expr(a), expr(b),expr(c),expr(d),expr(e),expr(f))) {
#  assign(x = paste(object, ".ref", sep = ""), value = eval(object))
# }
# save(a.ref,b.ref,c.ref,d.ref,e.ref,f.ref, file = "tests/testthat/ref_simulated_pheno_clines.Rda")
load("ref_simulated_pheno_clines.Rda")

test_that("clines generated by sim_pheno_cline are correct", {
  expect_equal(a, a.ref)
  expect_equal(b, b.ref)
  expect_equal(c, c.ref)
  expect_equal(d, d.ref)
  expect_equal(e, e.ref)
  expect_equal(f, f.ref)
})

