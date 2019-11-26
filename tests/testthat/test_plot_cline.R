context("plot_cline")

load("ref_geno_stanfit.Rda")
load("ref_geno_data.Rda")
load("ref_pheno_stanfit.Rda")
load("ref_simulated_geno_clines.Rda")
load("ref_pheno_data.Rda")


test_that("plot_cline checks args for type", {
  expect_error(plot_cline(stanfit = "ABC",
                               data = a.ref), "stanfit")
  expect_error(plot_cline(stanfit = ref_geno_stanfit,
                               data = "a.ref"), "data frame")
  expect_error(plot_cline(stanfit = ref_geno_stanfit,
                          data = a.ref, add.obs = "XXX"), "True")
  expect_error(plot_cline(stanfit = ref_geno_stanfit,
                          data = ref.geno.data, confidence = "XXX"), "TRUE")
  expect_error(plot_cline(stanfit = ref_geno_stanfit,
                          data = ref.geno.data, clear.cache = "XXX"), "TRUE")
  expect_error(plot_cline(ref_geno_stanfit, data = ref.geno.data, distance = 3, prob = "x"), "numeric")
  expect_error(plot_cline(stanfit = ref_geno_stanfit,
                               data = data.frame(transectDist = c("X", "Y"))), "numeric")
})

test_that("plot_cline checks args for appropriateness", {
  # Prob is good
  expect_error(plot_cline(ref_geno_stanfit, data = ref.geno.data, distance = 3, prob = c(0.6, 0.9)), "length")
  expect_error(plot_cline(ref_geno_stanfit, data = ref.geno.data, distance = 3, prob = -.5), "between")
  expect_error(plot_cline(ref_geno_stanfit, data = ref.geno.data, distance = 3, prob = 1.1), "between")
  # Colors are OK
  expect_error(plot_cline(stanfit = ref_geno_stanfit,
                          data = ref.geno.data, point.col = "xxxx"), "valid")
  expect_error(plot_cline(stanfit = ref_geno_stanfit,
                          data = ref.geno.data, cline.col = "xxxx"), "valid")
  # Check plotting extra args
  expect_error(plot_cline(ref_geno_stanfit, data = ref.geno.data, distance = 3, ann = "blue"), "default")
  expect_error(plot_cline(ref_geno_stanfit, data = ref.geno.data, distance = 3, border = "blue"), "default")
  expect_error(plot_cline(ref_geno_stanfit, data = ref.geno.data, distance = 3, col = "blue"), "default")
  expect_error(plot_cline(ref_geno_stanfit, data = ref.geno.data, distance = 3, type = "p"), "default")
  expect_error(plot_cline(ref_geno_stanfit, data = ref.geno.data, distance = 3, ylim = "blue"), "default")
})


test_that("plot_cline checks that necessary columns are present", {
  # For general plotting
  expect_error(plot_cline(stanfit = ref_geno_stanfit,
                          data = dplyr::select(ref.geno.data, -transectDist)), "contain")
  # When adding original data
  # geno
  expect_error(suppressWarnings(plot_cline(stanfit = ref_geno_stanfit,
                                           data = dplyr::select(ref.geno.data, -AA),
                                           add.obs = T)),
               "Necessary")
  # phno
  expect_error(suppressWarnings(plot_cline(stanfit = ref_pheno_stanfit,
                                           data = dplyr::select(ref.pheno.data, -traitValue),
                                           add.obs = T)),
               "Necessary")

})

test_that("plot_cline checks that stanfit and data frame match", {
  # geno
  expect_warning(plot_cline(stanfit = ref_pheno_stanfit,
                            data = ref.pheno.data[-1,]), "different")
  # pheno
  expect_warning(plot_cline(stanfit = ref_geno_stanfit,
                            data = b.ref), "different")

})

# # Output
# Don't think there's an easy way yet to check that the graphical output is correct
# Will want to check out the vdiffr package int he future.
# For now, can at least check that it makes it through and returns invisible NULL
test_that("plot_pheno_cline ouputs invisible NULL", {
  # geno
  expect_equal(suppressWarnings(plot_cline(stanfit = ref_geno_stanfit,
                                           data = ref.geno.data,
                                           add.obs = T)), NULL)
  expect_equal(suppressWarnings(plot_cline(stanfit = ref_geno_stanfit,
                                           data = ref.geno.data,
                                           add.obs = T,
                                           confidence = T,
                                           prob = 0.5)), NULL)
  # pheno
  expect_equal(suppressWarnings(plot_cline(stanfit = ref_pheno_stanfit,
                                           data = ref.pheno.data,
                                           add.obs = T)), NULL)
  expect_equal(suppressWarnings(plot_cline(stanfit = ref_pheno_stanfit,
                                           data = ref.pheno.data,
                                           add.obs = T,
                                           confidence = T,
                                           prob = 0.5)), NULL)
})

