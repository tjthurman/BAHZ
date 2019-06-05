context("prior checking")

test_that("assign_stan_dist_int returns proper values", {
  expect_equal(assign_stan_dist_int("normal"), 0)
  expect_equal(assign_stan_dist_int("uniform"), 1)
  expect_equal(assign_stan_dist_int("exponential"), 2)
  expect_equal(assign_stan_dist_int("beta"), 3)
  expect_error(assign_stan_dist_int("lkjhgkljhg"))
})

test_that("check_prior_specification works properly", {
  expect_true(check_prior_specification("normal", "normal(32,45)"))
  expect_true(check_prior_specification("uniform", "uniform(32,45)"))
  expect_true(check_prior_specification("exponential", "exponential(45)"))
  expect_false(check_prior_specification("normal", "normal(32)"))
  expect_false(check_prior_specification("uniform", "uniform(32)"))
  expect_false(check_prior_specification("exponential", "exponential(32,45)"))
  expect_false(check_prior_specification("beta", "beta(45)"))
})

test_that("check_prior_supported works properly", {
  # Supported
  for (param in c("center", "width")) {
    expect_true(check_prior_supported(param, "normal"))
    expect_true(check_prior_supported(param, "uniform"))
  }
  # pmin and pmax
  for (param in c("pmin", "pmax")) {
    expect_true(check_prior_supported(param, "uniform"))
    expect_true(check_prior_supported(param, "beta"))
    expect_true(check_prior_supported(param, "normal"))
  }
  # taus and f
  for (param in c("tauL", "tauR", "tauM", "f")) {
    expect_true(check_prior_supported(param, "uniform"))
    expect_true(check_prior_supported(param, "beta"))
  }
  #deltas
  for (param in c("deltaL", "deltaR", "deltaM")) {
    expect_true(check_prior_supported(param, "exponential"))
  }

  # Unsupported
  for (param in c("center", "width",
                  "pmin", "pmax",
                  "tauL", "tauR", "tauM", "f",
                  "deltaL", "deltaR", "deltaM")) {
    expect_false(check_prior_supported(param, "lkjhglkjhg"))
  }
})
