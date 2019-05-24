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
})
test_that("plot_geno_cline checks args for appropriateness", {
  expect_warning(plot_geno_cline(stanfit = ref_stanfit,
                                    data = a.ref), "different")
  expect_error(plot_geno_cline(stanfit = ref_stanfit,
                                    data = a.ref, point.col = "xxxx"), "valid")
})

# Output
# Don't think there's an easy way to check that the graphical output is correct
# but can atleast check that it makes it through and returns invisible NULL
test_that("plot_geno_cline ouputs invisible NULL", {
  expect_equal(suppressWarnings(plot_geno_cline(stanfit = ref_stanfit,
                                                data = a.ref,
                                                add.obs.freqs = T)), NULL)
})

