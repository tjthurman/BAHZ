context("init_single_chain")

# Argument checking
test_that("init_single_chain checks arguments correctly", {
  expect_error(init_single_chain("prior_config_test1.yaml", tails = "xxx"), "ind")
})

# Test output
test_that("init_single_chain outputs an expression", {
  expect_equal(class(init_single_chain("prior_config_test1.yaml", tails = "none")), "call")
})
test_that("init_single_chain outputes are of appropriate length for the different tails", {
  expect_equal(length(init_single_chain("prior_config_test1.yaml", tails = "left")), 7)
  expect_equal(length(init_single_chain("prior_config_test1.yaml", tails = "ind")), 9)
})
test_that("init_single_chain output is correct", {
  expect_equal_to_reference(init_single_chain("prior_config_test1.yaml", tails = "ind"),
                            file = "ref_init_single_chain_bi.Rda")
})
