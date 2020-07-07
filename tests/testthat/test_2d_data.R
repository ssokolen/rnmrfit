context("2D Basics")

#==============================================================================>
test_that("2d projection works", {

  nmroptions$sf <- 500

  p <- tibble(direct.shift = c(-1, 1), indirect.shift = c(-1, 1), 
              intensity = cmplx2(rr=c(-1, 1)))
  d2 <- new("NMRData2D")
  d2@acqus <- list(direct = list(sfo1 = nmroptions$sf),
                   indirect = list(sfo1 = nmroptions$sf))
  d2@procs <- list(direct = list(), indirect = list())
  d2@processed <- p

  p <- tibble(direct.shift = c(-1, 1), intensity = cmplx1(r=c(-1, 1)))
  d1 <- new("NMRData1D")
  d1@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d1@procs <- list(direct = list())
  d1@processed <- p

  # Checking getters
  expect_equal(direct(d2), d1)
})
