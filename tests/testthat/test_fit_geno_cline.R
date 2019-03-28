context("fit geno cline")

set.seed(18)
data <- sim_geno_cline(transect_distances = seq(-300,300,75), n_ind = 40, Fis = 0,
                       decrease = F, center = 10, width = 50, pmin = 0.03, pmax = .95)

# Test argument checking
test_that("fit geno cline checks arguments correctly", {
  expect_error(fit_geno_cline(data = data, prior_file = "prior_config_test1.yaml",
                              type = "xx", tails = "none", chains = 1), "multi")
  expect_error(fit_geno_cline(data = data, prior_file = "prior_config_test1.yaml",
                              type = "bi", tails = "xxx", chains = 1), "ind")
  expect_error(fit_geno_cline(data = data, prior_file = "prior_config_test1.yaml",
                              type = "bi", tails = "none", chains = "xxx"), "numeric")
})

# Test functionality
# Give it a really simple model to fit, and make sure the result is a stanfit.

z_p <- suppressWarnings(fit_geno_cline(data = data, prior_file = "prior_config_test1.yaml",
                      type = "bi", tails = "none", chains = 1))
test_that("fit geno cline runs the model and makes a stanfit object", {
  expect_equal(class(fit_geno_cline(data = data, prior_file = "prior_config_test1.yaml",
                                    type = "bi", tails = "none", chains = 1))[1], "stanfit")
})
z_p2 <- suppressWarnings(fit_geno_cline(data = data, prior_file = "prior_config_test1.yaml",
                                        type = "bi", tails = "none", chains = 1, control = list(adapt_delta = 0.85)))
test_that("fit geno cline allows stan control parameters to be changed", {
  expect_equal(attr(z_p2@sim$samples[[1]], "args")$control$adapt_delta, 0.85)
})

