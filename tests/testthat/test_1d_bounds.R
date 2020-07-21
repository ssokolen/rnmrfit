context("1D Bounds")

nmroptions$sf <- 500

################################################################################
# Resonance

object <- nmrresonance_1d("0.5 d 20", id = "R1")

#==============================================================================>
test_that("1d absolute general bounds work with resonance", {

  par <- list(position = c(2, 3), width = c(0.005, 5), 
              height = c(0, 10), fraction.gauss = c(0.1, 0.9), 
              widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par)
  
  check_general_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d widen general bounds work with resonance", {

  par.1 <- list(position = c(0, 0.5), width = c(0, 0.5), 
                height = c(0, 0.5), fraction.gauss = c(0.1, 0.9),
                widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par.1)

  par.2 <- list(position = c(-1, 1), width = c(-1, 1), 
                height = c(-1, 1), fraction.gauss = c(0, 1),
                widen = FALSE, object = object)

  object <- do.call(set_general_bounds, par.2)

  check_general_bounds(object, par.1, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d relative general bounds work with resonance", {

  d <- nmrdata_1d_from_scaffold(object)
  d@processed <- tibble(direct.shift = c(-1, 1), 
                        intensity = complex(re=c(-1, 1)))

  par <- list(position = c(0, 1), width = c(0, 1), 
              height = c(0, 1), nmrdata = d, widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par)

  par <- list(position = c(-1, 1), width = c(0, 2*nmroptions$sf), 
              height = c(0, 2))

  check_general_bounds(object, par, 1e-6)
})

#==============================================================================>
test_that("1d absolute offset bounds work with resonance", {

  par <- list(position = c(-0.05, 0.05), width = c(-0.5, 1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = FALSE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d relative offset bounds work with resonance", {

  par <- list(position = c(-0.01, 0.01), width = c(-0.1, 0.1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = TRUE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})



################################################################################
# Species

r1 <- nmrresonance_1d("0.3 d 20", id = "R1")
r2 <- nmrresonance_1d("0.7 d 20", id = "R2")
object <- nmrspecies_1d(list(r1, r2), id = "S1")

#==============================================================================>
test_that("1d absolute general bounds work with species", {

  par <- list(position = c(2, 3), width = c(0.005, 5), 
              height = c(0, 10), fraction.gauss = c(0.1, 0.9),
              widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par)

  check_general_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d widen general bounds work with species", {

  par.1 <- list(position = c(0, 0.5), width = c(0, 0.5), 
                height = c(0, 0.5), fraction.gauss = c(0.1, 0.9),
                widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par.1)

  par.2 <- list(position = c(-1, 1), width = c(-1, 1), 
                height = c(-1, 1), fraction.gauss = c(0, 1),
                widen = FALSE, object = object)

  object <- do.call(set_general_bounds, par.2)

  check_general_bounds(object, par.1, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d relative general bounds work with species", {

  d <- nmrdata_1d_from_scaffold(object)
  d@processed <- tibble(direct.shift = c(-1, 1), 
                        intensity = complex(re=c(-1, 1)))

  par <- list(position = c(0, 1), width = c(0, 1), 
              height = c(0, 1), nmrdata = d, widen = TRUE,
              object = object)

  object <- do.call(set_general_bounds, par)

  par <- list(position = c(-1, 1), width = c(0, 2*nmroptions$sf), 
              height = c(0, 2))

  check_general_bounds(object, par, 1e-6)
})

#==============================================================================>
test_that("1d absolute offset bounds work with species", {

  par <- list(position = c(-0.05, 0.05), width = c(-0.5, 1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = FALSE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d relative offset bounds work with species", {

  par <- list(position = c(-0.01, 0.01), width = c(-0.1, 0.1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = TRUE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})



################################################################################
# Fit

r1 <- nmrresonance_1d("0.3 d 20", id = "R1")
r2 <- nmrresonance_1d("0.7 d 20", id = "R2")
f <- nmrspecies_1d(list(r1, r2), id = "S1")

d <- nmrdata_1d_from_scaffold(f)
object <- nmrfit_1d(f, d, delay = TRUE)

#==============================================================================>
test_that("1d absolute general bounds work with fit", {

  par <- list(position = c(2, 3), width = c(0.005, 5), 
              height = c(0, 10), fraction.gauss = c(0.1, 0.9),
              widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par)

  check_general_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d widen general bounds work with fit", {

  par.1 <- list(position = c(0, 0.5), width = c(0, 0.5), 
                height = c(0, 0.5), fraction.gauss = c(0.1, 0.9),
                widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par.1)

  par.2 <- list(position = c(-1, 1), width = c(-1, 1), 
                height = c(-1, 1), fraction.gauss = c(0, 1),
                widen = FALSE, object = object)

  object <- do.call(set_general_bounds, par.2)

  check_general_bounds(object, par.1, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d relative general bounds work with fit", {

  d <- nmrdata_1d_from_scaffold(object)
  d@processed <- tibble(direct.shift = c(-1, 1), 
                        intensity = complex(re=c(-1, 1)))

  par <- list(position = c(0, 1), width = c(0, 1), 
              height = c(0, 1), nmrdata = d, widen = TRUE,
              object = object)

  object <- do.call(set_general_bounds, par)

  par <- list(position = c(-1, 1), width = c(0, 2*nmroptions$sf), 
              height = c(0, 2))

  check_general_bounds(object, par, 1e-6)
})

#==============================================================================>
test_that("1d absolute offset bounds work with fit", {

  par <- list(position = c(-0.05, 0.05), width = c(-0.5, 1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = FALSE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d relative offset bounds work with fit", {

  par <- list(position = c(-0.01, 0.01), width = c(-0.1, 0.1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = TRUE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})
