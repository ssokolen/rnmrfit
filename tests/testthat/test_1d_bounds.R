context("1D Bounds")

nmroptions$sf <- 500

################################################################################
# Resonance

#==============================================================================>
test_that("1d absolute general bounds work with resonance", {

  object <- gen_resonance_1d()

  par <- list(position = c(2, 3), width = c(0.005, 5), 
              height = c(0, 10), fraction.gauss = c(0.1, 0.9), 
              widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par)
  
  check_general_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d widen general bounds work with resonance", {

  object <- gen_resonance_1d()

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

  object <- gen_resonance_1d()

  p <- tibble(direct.shift = c(-1, 1), intensity = complex(re=c(-1, 1)))
  d <- new("NMRData1D")
  d@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  par <- list(position = c(0, 1), width = c(0, 1), 
              height = c(0, 1), nmrdata = d, widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par)

  par <- list(position = c(-1, 1), width = c(0, 2*nmroptions$sf), 
              height = c(0, 2))

  check_general_bounds(object, par, 1e-6)
})

#==============================================================================>
test_that("1d absolute offset bounds work with resonance", {

  nmroptions$sf <- 500

  object <- gen_resonance_1d()

  par <- list(position = c(-0.05, 0.05), width = c(-0.5, 1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = FALSE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d relative offset bounds work with resonance", {

  nmroptions$sf <- 500

  object <- gen_resonance_1d()

  par <- list(position = c(-0.01, 0.01), width = c(-0.1, 0.1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = TRUE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})



################################################################################
# Species



#==============================================================================>
test_that("1d absolute general bounds work with species", {

  object <- gen_species_1d()

  par <- list(position = c(2, 3), width = c(0.005, 5), 
              height = c(0, 10), fraction.gauss = c(0.1, 0.9),
              widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par)

  check_general_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d widen general bounds work with species", {

  object <- gen_species_1d()

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

  object <- gen_species_1d()

  p <- tibble(direct.shift = c(-1, 1), intensity = complex(re=c(-1, 1)))
  d <- new("NMRData1D")
  d@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

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

  object <- gen_species_1d()

  par <- list(position = c(-0.05, 0.05), width = c(-0.5, 1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = FALSE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d relative offset bounds work with species", {

  object <- gen_species_1d()

  par <- list(position = c(-0.01, 0.01), width = c(-0.1, 0.1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = TRUE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})



################################################################################
# Fit



#==============================================================================>
test_that("1d absolute general bounds work with fit", {

  object <- gen_species_1d() %>% gen_fit_1d()

  par <- list(position = c(2, 3), width = c(0.005, 5), 
              height = c(0, 10), fraction.gauss = c(0.1, 0.9),
              widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par)

  check_general_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d widen general bounds work with fit", {

  object <- gen_species_1d() %>% gen_fit_1d()

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

  object <- gen_species_1d() %>% gen_fit_1d()

  p <- tibble(direct.shift = c(-1, 1), intensity = complex(re=c(-1, 1)))
  d <- new("NMRData1D")
  d@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

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

  object <- gen_species_1d() %>% gen_fit_1d()

  par <- list(position = c(-0.05, 0.05), width = c(-0.5, 1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = FALSE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})

#------------------------------------------------------------------------------
test_that("1d relative offset bounds work with fit", {

  object <- gen_species_1d() %>% gen_fit_1d()

  par <- list(position = c(-0.01, 0.01), width = c(-0.1, 0.1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = TRUE,
              object = object)

  object <- do.call(set_offset_bounds, par)

  check_offset_bounds(object, par, 1e-6)
})
