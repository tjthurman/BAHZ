context("prep_geno_data")

# Test argument checking
test_that("prep_geno_data checks arguments", {
  expect_error(prep_geno_data(as.matrix(c(1,2,3)), type = "bi"), "tibble")
  expect_error(prep_geno_data(as.list(c(1,2,3)), type = "bi"), "tibble")
  expect_error(prep_geno_data(data.frame(nFocalAllele = c(1,2,3),
                                         nTotalAlleles = c(1,2,3),
                                         transectDist = c(1,2,3)), type = "hxy"), "either")
})

test_that("prep_geno_data checks for proper columns", {
  expect_error(prep_geno_data(data.frame(nFocalAllele = c(1,2,3),
                                         nTotalAlleles = c(1,2,3),
                                         transectDist = c(1,2,3)), type = "multi"), "multinomial")
  expect_error(prep_geno_data(data.frame(nFocalAllele = c(1,2,3),
                                         nTotalalleles = c(1,2,3),
                                         transectDist = c(1,2,3)), type = "bi"), "binomial")
  expect_error(prep_geno_data(data.frame(nFocalAllele = c(1,2,3),
                                         nTotalAlleles = c(1,2,3),
                                         transectDist = c("1","2","3")), type = "bi"), "numeric")
  expect_error(prep_geno_data(data.frame(nFocalAllele = c(1,2,3),
                                         nTotalAlleles = c(1,2.5,3),
                                         transectDist = c(1,2,3)), type = "bi"), "integer")
  expect_error(prep_geno_data(data.frame(AA = c(1,2,3),
                                         Aa = c(1,2.5,3),
                                         aa = c(1,2.5,3),
                                         transectDist = c(1,2,3)), type = "multi"), "integer")
})

test_that("prep_geno_data handles flat cline error", {
  expect_error(prep_geno_data(data.frame(nFocalAllele = as.integer(c(10,5,10)),
                                         nTotalAlleles = as.integer(c(14,10,14)),
                                         transectDist = c(1,2,3)), type = "bi"), "increasing")
  expect_error(prep_geno_data(data.frame(AA = as.integer(c(10,5,10)),
                                         Aa = as.integer(c(1,10,1)),
                                         aa = as.integer(c(0,10,0)),
                                         transectDist = c(1,2,3)), type = "bi"), "increasing")
  expect_error(prep_geno_data(data.frame(AA = as.integer(c(10,5,10)),
                                         Aa = as.integer(c(1,10,1)),
                                         aa = as.integer(c(0,10,0)),
                                         transectDist = c(1,2,3)), type = "multi"), "increasing")
})

# Test output
test_that("prep_geno_data gives correct output", {
   expect_equal_to_reference(prep_geno_data(data.frame(nFocalAllele = as.integer(c(10,5,2)),
                                                       nTotalAlleles = as.integer(c(14,10,14)),
                                                       transectDist = c(1,2,3)), type = "bi"), file = "ref_prep_geno_data1.Rda")
  expect_equal_to_reference(prep_geno_data(data.frame(AA = as.integer(c(10,5,2)),
                                                      Aa = as.integer(c(1,10,1)),
                                                      aa = as.integer(c(0,10,0)),
                                                      transectDist = c(1,2,3)), type = "bi"), file = "ref_prep_geno_data2.Rda")
  expect_equal_to_reference(prep_geno_data(data.frame(AA = as.integer(c(10,5,2)),
                                                      Aa = as.integer(c(1,10,1)),
                                                      aa = as.integer(c(0,10,0)),
                                                      transectDist = c(1,2,3)), type = "multi"), file = "ref_prep_geno_data3.Rda")
})




