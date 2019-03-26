context("parse prior file")

test_that("parse prior checks number, order, and names of priors", {
  expect_error(parse_prior_file("prior_config_wrong_number.yaml"), "number")
  expect_error(parse_prior_file("prior_config_wrong_order.yaml"), "invalid")
  expect_error(parse_prior_file("prior_config_wrong_name.yaml"), "xyz")
})

test_that("parse prior reads in file properly", {
  expect_equal_to_reference(parse_prior_file("prior_config_test1.yaml"), file = "ref_parse_prior.Rda")
})
