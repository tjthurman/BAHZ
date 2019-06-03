context("prep_pheno_data")

# Test argument checking
test_that("prep_pheno_data checks arguments", {
  expect_error(prep_pheno_data(as.matrix(c(1,2,3))), "tibble")
  expect_error(prep_pheno_data(as.list(c(1,2,3))), "tibble")
})

test_that("prep_pheno_data checks for proper columns", {
  expect_error(prep_pheno_data(data.frame(traitValue = c(1,2,3),
                                          transectDist = c("A",2,3)), "numeric"))
  expect_error(prep_pheno_data(data.frame(traitValue = c("A",2,3),
                                          transectDist = c(1,2,3)), "numeric"))
  expect_error(prep_pheno_data(data.frame(traitValue = c(1,2,3),
                                          transectdist = c(1,2,3)), "transectDist"))

})

test_that("prep_pheno_data handles flat cline error", {
  expect_error(prep_pheno_data(data.frame(traitValue = c(1,2,1),
                                          transectDist = c(1,2,3)), "equal"))
})

# Test output
test_that("prep_pheno_data gives correct output", {
  expect_equal_to_reference(prep_pheno_data(data.frame(traitValue = c(10,15,20),
                                                       transectDist = c(1,2,3))),
                                            file = "ref_prep_pheno_data1.Rda")
  expect_equal_to_reference(prep_pheno_data(data.frame(traitValue = c(15,20,10),
                                                       transectDist = c(2,3,1))),
                            file = "ref_prep_pheno_data1.Rda")
})




