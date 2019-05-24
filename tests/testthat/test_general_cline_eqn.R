context("general cline equation")

# Test argument checking
test_that("general cline equation checks arguments for type", {
  expect_error(general_cline_eqn(transectDist = as.matrix(c(1,2,3))), "vector")
  expect_error(general_cline_eqn(transectDist = "hello"), "numeric")
  expect_error(general_cline_eqn(transectDist = 20, decrease = "xxx"), "decrease")
  expect_error(general_cline_eqn(transectDist = 20, decrease = T,
                                 center = as.matrix(c(1,2,3)), width = 10), "vector")
  expect_error(general_cline_eqn(transectDist = 20, decrease = T,
                                 center = 20, width = "10"), "numeric")
  expect_error(general_cline_eqn(transectDist = 20, decrease = T,
                                 center = 20, width = 10,
                                 pmin = 0.2, pmax = 0.8, deltaL = c(1,2)), "length")
})

test_that("general cline equation checks arguments for appropriateness", {
  expect_error(general_cline_eqn(transectDist = 20, decrease = T,
                                 center = 20, width = -10,
                                 pmin = 0.2, pmax = 0.8), "width")
  expect_error(general_cline_eqn(transectDist = 20, decrease = T,
                                 center = 20, width = 10,
                                 pmin = -0.2, pmax = 0.8), "between")
  expect_error(general_cline_eqn(transectDist = 20, decrease = T,
                                 center = 20, width = 10,
                                 pmin = 0.2, pmax = 0.8, tauL = 1.1, deltaL = 10), "between")
  expect_error(general_cline_eqn(transectDist = 20, decrease = T,
                                 center = 20, width = 10,
                                 pmin = 0.2, pmax = 0.8, tauL = 1), "tauL")
  expect_error(general_cline_eqn(transectDist = 20, decrease = T,
                                 center = 20, width = 10,
                                 pmin = 0.2, pmax = 0.8, deltaR = 1), "deltaR")

})

# Test results
test_that("general_cline_equation gives the proper results", {
  expect_equal(general_cline_eqn(transectDist = 70, decrease = F,
                                 center = 20, width = 10,
                                 pmin = 0.2, pmax = 0.8), 0.8)
  expect_equal(general_cline_eqn(transectDist = 70, decrease = T,
                                 center = 20, width = 10,
                                 pmin = 0.2, pmax = 0.8), 1-0.8)
  expect_equal(general_cline_eqn(transectDist = 20, decrease = F,
                                 center = 20, width = 10), 0.5)
  expect_equal(general_cline_eqn(transectDist = 85, decrease = T,
                                 center = 100, width = 60,
                                 pmin = 0.2, pmax = 0.8, deltaL = 10, tauL = 0.6), 1- 0.3783489, tolerance = 1e-7)
  expect_equal(general_cline_eqn(transectDist = 115, decrease = F,
                                 center = 100, width = 60,
                                 pmin = 0.2, pmax = 0.8, deltaR = 10, tauR = 0.6), 1- 0.3783489, tolerance = 1e-7)
  expect_equal(general_cline_eqn(transectDist = 115, decrease = F,
                                 center = 100, width = 60,
                                 pmin = 0.2, pmax = 0.8,
                                 deltaR = 10, tauR = 0.6,
                                 deltaL = 14, tauL = 0.5), 0.6216511, tolerance = 1e-7)
  expect_equal(general_cline_eqn(transectDist = 85, decrease = F,
                                 center = 100, width = 60,
                                 pmin = 0.2, pmax = 0.8,
                                 deltaR = 10, tauR = 0.6,
                                 deltaL = 14, tauL = 0.5), 0.3653458, tolerance = 1e-7)
  expect_equal(general_cline_eqn(transectDist = c(85,115), decrease = F,
                                 center = 100, width = 60,
                                 pmin = 0.2, pmax = 0.8,
                                 deltaR = 10, tauR = 0.6,
                                 deltaL = 14, tauL = 0.5), c(0.3653458, 0.6216511), tolerance = 1e-7)

})

test_that("general_cline_eqn gives proper result at edge case where exp goes to infinity", {
  expect_equal(general_cline_eqn(transectDist = 120, decrease = F,
                                 center = 0, width = 2, deltaR = 1, tauR = 0.5), 1)

})

