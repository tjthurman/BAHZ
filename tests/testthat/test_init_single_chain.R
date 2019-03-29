context("init_single_chain")

# Argument checking
test_that("init_single_chain checks arguments correctly", {
  expect_error(init_single_chain("prior_config_test1.yaml", tails = "xxx", type = "bi"), "ind")
  expect_error(init_single_chain("prior_config_test1.yaml", tails = "none", type = "xxx"), "bi")
})

# Test output
test_that("init_single_chain outputs an expression", {
  expect_equal(class(init_single_chain("prior_config_test1.yaml", tails = "none", type = "bi")), "call")
})
test_that("init_single_chain outputes are of appropriate length for the different tails", {
  expect_equal(length(init_single_chain("prior_config_test1.yaml", tails = "left", type = "bi")), 7)
  expect_equal(length(init_single_chain("prior_config_test1.yaml", tails = "ind", type = "bi")), 9)
  expect_equal(length(init_single_chain("prior_config_test1.yaml", tails = "left", type = "multi")), 8)
  expect_equal(length(init_single_chain("prior_config_test1.yaml", tails = "ind", type = "multi")), 10)
})
test_that("init_single_chain output is correct", {
  expect_equal_to_reference(init_single_chain("prior_config_test1.yaml", tails = "ind", type = "bi"),
                            file = "ref_init_single_chain_bi.Rda")
  expect_equal_to_reference(init_single_chain("prior_config_test1.yaml", tails = "ind", type = "multi"),
                            file = "ref_init_single_chain_multi.Rda")
})
