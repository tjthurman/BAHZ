context("cline summary")

load("ref_stanfit.Rda")

# Test of argument checks
test_that("cline summary argument checks work", {
  expect_error(cline_summary(stanfit = "a"), "stanfit")
  expect_error(cline_summary(ref_stanfit, prob = "x"), "numeric")
  expect_error(cline_summary(ref_stanfit, prob = c(0.6, 0.9)), "length")
  expect_error(cline_summary(ref_stanfit, prob = -.5), "between")
  expect_error(cline_summary(ref_stanfit, prob = 1.1), "between")
  expect_error(cline_summary(ref_stanfit, method = "xxx"), "method")
  expect_error(cline_summary(ref_stanfit, show.all = "xxx"), "show.all")
})


# Test of functionality
test_that("cline summary provides correct summary", {
  expect_equal_to_reference(cline_summary(ref_stanfit,
                                          prob = 0.5), file = "ref_summary1.Rda")
  expect_equal_to_reference(cline_summary(ref_stanfit, method = "ET",
                                          prob = 0.875), file = "ref_summary2.Rda")
  expect_equal_to_reference(cline_summary(ref_stanfit, show.all = T,
                                          prob = 0.65), file = "ref_summary3.Rda")
})
