context("make_prior_config")

# Test argument checking
test_that("make prior config checks arguments", {
  expect_error(make_prior_config(name = c("one", "two")), "length")
  expect_error(make_prior_config(name = 2), "character")
  expect_error(make_prior_config(name = "prior_config.xyz"), "yaml")
  expect_error(make_prior_config(name = "prior_config.yaml", overwrite = "X"), "overwrite")
})

make_prior_config(name = "prior_config_template_test.yaml")
# Test functionality
test_that("make prior config makes the file", {
  expect_true(file.exists("prior_config_template_test.yaml"))
})
test_that("make prior config won't overwrite by default", {
  expect_error(make_prior_config(name = "prior_config_template_test.yaml"), "overwrite")
})
file.remove("prior_config_template_test.yaml")
# run it, check that file exists in test directory?
