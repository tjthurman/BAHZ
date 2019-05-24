context("predict_geno_cline")

load("ref_stanfit.Rda")
load("ref_simulated_geno_clines.Rda")

#Argument checking
test_that("predict_geno_cline checks args for type", {
  expect_error(predict_geno_cline(stanfit = "ABC",
                                  distance = 3), "stanfit")
  expect_error(predict_geno_cline(stanfit = ref_stanfit,
                                  distance = as.matrix(3)), "vector")
  expect_error(predict_geno_cline(stanfit = ref_stanfit,
                                  distance = "X"), "numeric")
  })

# Output checking
# Right now only checks proper output for
# Simplest model (no tails)
test_that("predict_geno_cline outputs data correctly", {
  expect_true(is.data.frame(predict_geno_cline(stanfit = ref_stanfit,
                                               distance = 3)))
  expect_equal_to_reference(predict_geno_cline(stanfit = ref_stanfit,
                                               distance = 0:100), file = "ref_pred_geno_cline1.Rda")
  expect_true(dim(predict_geno_cline(stanfit = ref_stanfit,
                                     distance = seq(0, 10, length.out = 7)))[1] == 7)
})
