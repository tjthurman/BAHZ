context("prep prior list")

test_that("prep pior list gives informative warnings when priors are wrong", {
  expect_error(prep_prior_list("prior_config_test_unsupported_dist.yaml"), "supported")
  expect_error(prep_prior_list("prior_config_test_misspecified_dist.yaml"), "correct")
})

test_that("prep prior list makes the right list", {
  expect_equal_to_reference(prep_prior_list("prior_config_test1.yaml"), file = "ref_prep_prior_list1.Rda")
  expect_equal_to_reference(prep_prior_list("prior_config_test2.yaml"), file = "ref_prep_prior_list2.Rda")
})
