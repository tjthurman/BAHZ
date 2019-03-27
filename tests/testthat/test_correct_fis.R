context("correct_fis")

error.column <- data.frame(a = "1", b = "c")
test_that("correct fis throw error without fis column", {
  expect_error(correct_fis(error.column), "named")
})

na.handle <- data.frame(Fis = c(NaN, NaN, NaN))
neg.handle <- data.frame(Fis = c(-19, -0.00048, -1))
result <- data.frame(Fis = c(0, 0, 0))

test_that("correct fis handles NaNs", {
  expect_equal(correct_fis(na.handle), result)
})
test_that("correct fis handles negatives", {
  expect_equal(correct_fis(neg.handle), result)
})
