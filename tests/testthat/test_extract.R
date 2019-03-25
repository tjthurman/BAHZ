context("extract functions")

test_that("extract functions are situation specific", {
  #First
  expect_error(extract_first("(8)"))
  #Last
  expect_error(extract_last("(8)"))
  #Error
  expect_error(extract_only("(8,8)"), "comma")
})


test_that("extract functions handle decimals", {
  # First
  expect_equal(extract_first("(0.08, 0.2)"), 0.08)
  expect_error(extract_first("(0..08, 0.2)"), "decimal")
  expect_equal(extract_first("(8, 0.2)"), 8)
  # Last
  expect_equal(extract_last("(0.08, 0.2)"), 0.2)
  expect_error(extract_last("(0.08, 0..2)"), "decimal")
  expect_equal(extract_last("(8, 2)"), 2)
  # Only
  expect_equal(extract_only("(0.08)"), 0.08)
  expect_error(extract_only("(0..08)"),  "decimal")
  expect_equal(extract_only("(8)"), 8)
})

test_that("extract functions handle non-numerics", {
  expect_error(extract_first("(XX, 0.01)"), "coerce")
  expect_error(extract_last("(1, XX)"), "coerce")
  expect_error(extract_only("(XX)"), "coerce")
})

test_that("extract functions handle whitespace", {
  # First
  expect_equal(extract_first("   (   08 , 2)"), 8)
  expect_equal(extract_first("   (   0.8 , 2)"), 0.8)
  expect_equal(extract_first("   (   0 8 , 2)"), 8)
  # Last
  expect_equal(extract_last("   (   08 ,   2   )"), 2)
  expect_equal(extract_last("   (   0.8 , .2)"), 0.2)
  expect_equal(extract_last("   (   0 8 , 2 0 )"), 20)
  # Only
  expect_equal(extract_only("   (   2     )"), 2)
  expect_equal(extract_only("   (   .2    )"), 0.2)
  expect_equal(extract_only("   ( 2 0 )"), 20)
})

test_that("extract functions handle negative values", {
  #First
  expect_equal(extract_first("(-0.03, 10)"), -0.03)
  expect_equal(extract_first("( - 3 , 10)"), -3)
  #Last
  expect_equal(extract_last("(-0.03, -10)"), -10)
  expect_equal(extract_last("( - 3 , - 1 0 . 0 0)"), -10)
})

