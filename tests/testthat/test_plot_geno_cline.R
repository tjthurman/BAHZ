context("plot_geno_cline")

load("ref_stanfit.Rda")
load("ref_simulated_geno_clines.Rda")

#Argument checking
test_that("plot_geno_cline checks args for type", {
  expect_error(plot_geno_cline(stanfit = "ABC",
                                  data = a.ref), "stanfit")
  expect_error(plot_geno_cline(stanfit = ref_stanfit,
                                  data = "a.ref"), "data frame")
  expect_error(plot_geno_cline(stanfit = ref_stanfit,
                                  data = dplyr::select(a.ref, -transectDist)), "contain")
  expect_error(suppressWarnings(plot_geno_cline(stanfit = ref_stanfit,
                               data = dplyr::select(a.ref, -AA),
                               add.obs.freqs = T)), "Necessary")
  expect_error(plot_geno_cline(stanfit = ref_stanfit,
                                  data = data.frame(transectDist = c("X", "Y"))), "numeric")
  expect_error(plot_geno_cline(stanfit = ref_stanfit,
                               data = a.ref, add.obs.freq = "XXX"), "True")
  expect_error(plot_geno_cline(stanfit = ref_stanfit,
                               data = a.ref, confidence = "XXX"), "TRUE")
  expect_error(plot_geno_cline(ref_stanfit, data = a.ref, distance = 3, prob = "x"), "numeric")
  expect_error(plot_geno_cline(ref_stanfit, data = a.ref, distance = 3, col = "blue"), "default")
  expect_error(plot_geno_cline(ref_stanfit, data = a.ref, distance = 3, type = "p"), "default")
})
test_that("plot_geno_cline checks args for appropriateness", {
  expect_warning(plot_geno_cline(stanfit = ref_stanfit,
                                    data = a.ref), "different")
  expect_error(plot_geno_cline(stanfit = ref_stanfit,
                                    data = a.ref, point.col = "xxxx"), "valid")
  expect_error(plot_geno_cline(stanfit = ref_stanfit,
                               data = a.ref, cline.col = "xxxx"), "valid")
  expect_error(plot_geno_cline(ref_stanfit, data = a.ref, distance = 3, prob = c(0.6, 0.9)), "length")
  expect_error(plot_geno_cline(ref_stanfit, data = a.ref, distance = 3, prob = -.5), "between")
  expect_error(plot_geno_cline(ref_stanfit, data = a.ref, distance = 3, prob = 1.1), "between")
})

# Output
# Don't think there's an easy way to check that the graphical output is correct
# but can atleast check that it makes it through and returns invisible NULL
test_that("plot_geno_cline ouputs invisible NULL", {
  expect_equal(suppressWarnings(plot_geno_cline(stanfit = ref_stanfit,
                                                data = a.ref,
                                                add.obs.freqs = T)), NULL)
  expect_equal(suppressWarnings(plot_geno_cline(stanfit = ref_stanfit,
                                                data = a.ref,
                                                add.obs.freqs = T,
                                                confidence = T,
                                                prob = 0.5)), NULL)
})

