context("predict_geno_cline")

load("ref_stanfit.Rda")
load("ref_simulated_geno_clines.Rda")

#Argument checking
test_that("predict_geno_cline checks args for type", {
  expect_error(predict_geno_cline(stanfit = "ABC",
                                  data = a.ref), "stanfit")
  expect_error(predict_geno_cline(stanfit = ref_stanfit,
                                  data = "a.ref"), "data frame")
  expect_error(predict_geno_cline(stanfit = ref_stanfit,
                                  data = dplyr::select(a.ref, -transectDist)), "contain")
  expect_error(predict_geno_cline(stanfit = ref_stanfit,
                                  data = data.frame(transectDist = c("X", "Y"))), "numeric")
  expect_error(predict_geno_cline(stanfit = ref_stanfit,
                                    data = a.ref, num.out = "XXX"), "numeric")
  })
test_that("predict_geno_cline checks args for appropriateness", {
  expect_warning(predict_geno_cline(stanfit = ref_stanfit,
                                    data = a.ref), "different")
})

# Output checking
# Right now only checks proper output for
# Simplest model (no tails)
test_that("predict_geno_cline outputs data correctly", {
  expect_true(is.data.frame(suppressWarnings(predict_geno_cline(stanfit = ref_stanfit,
                                                  data = a.ref))))
  expect_equal_to_reference(suppressWarnings(predict_geno_cline(stanfit = ref_stanfit,
                                               data = a.ref)), file = "ref_pred_geno_cline1.Rda")
  expect_true(dim(suppressWarnings(predict_geno_cline(stanfit = ref_stanfit,
                                 data = a.ref, num.out = 78)))[1] == 78)
  expect_equal_to_reference(suppressWarnings(predict_geno_cline(stanfit = ref_stanfit,
                                                                data = a.ref)), file = "ref_pred_geno_cline2.Rda")
})
