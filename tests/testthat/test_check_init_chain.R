context("check_init_chain")

width.bad <- list(center = 10,
                  width = -10,
                  pmin = 0.2,
                  pmax = 0.8,
                  deltaL = 2,
                  tauL = 0.5,
                  deltaR = 3,
                  tauR = 0.5)
pmin.bad <- list(center = 10,
                  width = 10,
                  pmin = -0.2,
                  pmax = 0.8,
                  deltaL = 2,
                  tauL = 0.5,
                  deltaR = 3,
                  tauR = 0.5)
pmax.bad <- list(center = 10,
                 width = 10,
                 pmin = 0.2,
                 pmax = 1.1,
                 deltaL = 2,
                 tauL = 0.5,
                 deltaR = 3,
                 tauR = 0.5)
pmax.small <- list(center = 10,
                    width = 10,
                    pmin = 0.2,
                    pmax = 0.1,
                    deltaL = 2,
                    tauL = 0.5,
                    deltaR = 3,
                    tauR = 0.5)
deltaL.bad <- list(center = 10,
                   width = 10,
                   pmin = 0.2,
                   pmax = 0.3,
                   deltaL = -2,
                   tauL = 0.5,
                   deltaR = 3,
                   tauR = 0.5)
deltaR.bad <- list(center = 10,
                   width = 10,
                   pmin = 0.2,
                   pmax = 0.3,
                   deltaL = 2,
                   tauL = 0.5,
                   deltaR = -3,
                   tauR = 0.5)
deltaM.bad <- list(center = 10,
                   width = 10,
                   pmin = 0.2,
                   pmax = 0.3,
                   deltaM = -2,
                   tauM = 0.5)
tauL.bad <- list(center = 10,
                 width = 10,
                 pmin = 0.2,
                 pmax = 0.3,
                 deltaL = 2,
                 tauL = -0.5,
                 deltaR = 3,
                 tauR = 0.5)
tauR.bad <- list(center = 10,
                 width = 10,
                 pmin = 0.2,
                 pmax = 0.3,
                 deltaL = 2,
                 tauL = 0.5,
                 deltaR = 3,
                 tauR = 1.1)
tauM.bad <- list(center = 10,
                 width = 10,
                 pmin = 0.2,
                 pmax = 0.3,
                 deltaL = 2,
                 tauM = -0.5,
                 deltaR = 3,
                 tauR = 0.5)
f.bad <- list(center = 10,
              width = 10,
              pmin = 0.2,
              pmax = 0.3,
              deltaL = 2,
              tauL = 0.5,
              deltaR = 3,
              tauR = 0.5,
              f = -1)
good <- list(center = 10,
              width = 10,
              pmin = 0.2,
              pmax = 0.3,
              deltaL = 2,
              tauL = 0.5,
              deltaR = 3,
              tauR = 0.5,
              f = 0.5)

test_that("check_init_chain output is correct", {
  expect_equal(check_init_chain(width.bad), "width")
  expect_equal(check_init_chain(pmin.bad), "pmin")
  expect_equal(check_init_chain(pmax.bad), "pmax")
  expect_equal(check_init_chain(pmax.small), c("pmin", "pmax"))
  expect_equal(check_init_chain(deltaL.bad), "deltaL")
  expect_equal(check_init_chain(deltaR.bad), "deltaR")
  expect_equal(check_init_chain(deltaM.bad), "deltaM")
  expect_equal(check_init_chain(tauL.bad), "tauL")
  expect_equal(check_init_chain(tauR.bad), "tauR")
  expect_equal(check_init_chain(tauM.bad), "tauM")
  expect_equal(check_init_chain(f.bad), "f")
  expect_null(check_init_chain(good))
})
