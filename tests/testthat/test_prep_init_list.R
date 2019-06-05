context("prep_init_list")

# Test argument checking
test_that("prep_init_list checks arguments", {
  expect_error(prep_init_list("prior_config_test1.yaml", tails = "xxx", chains = as.integer(3), type = "geno"), "ind")
  expect_error(prep_init_list("prior_config_test1.yaml", tails = "none", chains = "a", type = "geno"), "integer")
  expect_error(prep_init_list("prior_config_test1.yaml", tails = "none", chains = as.integer(c(1,2)), type = "geno"), "single")
  expect_error(prep_init_list("prior_config_test1.yaml", tails = "none", chains = as.integer(3), type = "xxx"), "should")
})

# Test functionality
test_that("prep_init_list fails when it can't initialize due to bad values", {
  expect_error(prep_init_list("prior_config_test_no_initial.yaml", tails = "none", chains = as.integer(1), type = "geno"), "Could not generate")
})
test_that("prep_init_list returns the proper values", {
  set.seed(763)
  expect_equal_to_reference(prep_init_list("prior_config_test1.yaml",
                                           tails = "ind",
                                           chains = as.integer(5),
                                           type = "pheno"), file = "ref_prep_init_list_bi.Rda")
})
