context("fit pheno cline")

set.seed(22)
data <- sim_pheno_cline(transect_distances = seq(-300,300,75), n_ind = 30, sigma = 6,
                       decrease = F, center = 10, width = 35, pmin = 30, pmax = 50)

# Test argument checking for variance arg
test_that("fit pheno cline checks arguments", {
  expect_error(fit_pheno_cline(data = data,
                               prior_file = "prior_config_test3.yaml",
                               chains = 1, pheno_variance = "XXX"), "constant")
  expect_error(fit_pheno_cline(data = data, prior_file = "prior_config_test3.yaml",
                               chains = "xxx"), "numeric")
})


# Test functionality
# Give it a really simple model to fit, and make sure the result is a stanfit.

# z_p <- suppressWarnings(fit_pheno_cline(data = data, prior_file = "prior_config_test1.yaml",
#                       type = "bi", tails = "none", chains = 1))
test_that("fit pheno cline runs the model and makes a stanfit object", {
  expect_equal(class(suppressWarnings(fit_pheno_cline(data = data, prior_file = "prior_config_test3.yaml"
                                                      , chains = 1)))[1], "stanfit")
  expect_equal(class(suppressWarnings(fit_pheno_cline(data = data, prior_file = "prior_config_test3.yaml"
                                                      , chains = 1, pheno_variance = "pooled")))[1], "stanfit")
})

test_that("fit_pheno_cline works with user-defined init list", {
  expect_equal(class(suppressWarnings(fit_pheno_cline(data = data, prior_file = "prior_config_test3.yaml",
                                                     chains = 1,
                                                     init = list(list(center = 9,
                                                                 width =  49,
                                                                 pmin = 29,
                                                                 pmax = 49)))))[1], "stanfit")
})

z_p2 <- suppressWarnings(fit_pheno_cline(data = data, prior_file = "prior_config_test1.yaml",
                                        chains = 1, control = list(adapt_delta = 0.85)))
test_that("fit pheno cline allows stan control parameters to be changed", {
  expect_equal(attr(z_p2@sim$samples[[1]], "args")$control$adapt_delta, 0.85)
})
