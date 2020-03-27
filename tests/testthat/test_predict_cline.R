context("predict_cline")

load("ref_geno_stanfit.Rda")
load("ref_simulated_geno_clines.Rda")

#Argument checking
test_that("predict_cline checks args for type", {
  expect_error(predict_cline(stanfit = "ABC",
                                  distance = 3), "stanfit")
  expect_error(predict_cline(stanfit = ref_geno_stanfit,
                                  distance = as.matrix(3)), "vector")
  expect_error(predict_cline(stanfit = ref_geno_stanfit,
                                  distance = "X"), "numeric")

  expect_error(predict_cline(ref_geno_stanfit, distance = 3, prob = "x"), "numeric")
  expect_error(predict_cline(ref_geno_stanfit, distance = 3, prob = c(0.6, 0.9)), "length")
  expect_error(predict_cline(ref_geno_stanfit, distance = 3, prob = -.5), "between")
  expect_error(predict_cline(ref_geno_stanfit, distance = 3, prob = 1.1), "between")
  expect_error(predict_cline(ref_geno_stanfit, distance = 3, confidence = "XXX"), "TRUE")
  expect_error(predict_cline(ref_geno_stanfit, distance = 3, progress = "XXX"), "TRUE")
  expect_error(predict_cline(ref_geno_stanfit, distance = 3, method = "XXX"), "HPDI")
  })

# Output checking
# Right now only checks proper output for
# Simplest model (no tails)
test_that("predict_cline outputs data correctly", {
  expect_true(is.data.frame(predict_cline(stanfit = ref_geno_stanfit,
                                               distance = 3)))
  expect_true(is.data.frame(predict_cline(stanfit = ref_geno_stanfit,
                                          distance = 3, confidence = T)))
  expect_equal(dim(predict_cline(stanfit = ref_geno_stanfit,
                                          distance = 3, confidence = T))[2], 5)
  expect_equal_to_reference(predict_cline(stanfit = ref_geno_stanfit,
                                               distance = 0:100),
                            file = "ref_pred_cline1.Rda")
  expect_equal_to_reference(predict_cline(stanfit = ref_geno_stanfit,
                                          distance = 0:100),
                            confidence = T,
                            file = "ref_pred_cline2.Rda")
  expect_equal_to_reference(predict_cline(stanfit = ref_geno_stanfit,
                                          distance = 0:100,
                                          method = "ET"),
                            confidence = T,
                            file = "ref_pred_cline3.Rda")
  expect_true(dim(predict_cline(stanfit = ref_geno_stanfit,
                                     distance = seq(0, 10, length.out = 7)))[1] == 7)
})
