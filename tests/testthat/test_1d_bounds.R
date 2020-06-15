context("1D Bounds")

nmroptions$sf <- 500

#==============================================================================>
# The same tests wll be run on resonance/species/fit
check_general_bounds <- function(object, par_in, par_out, tolerance) {

  testthat_tolerance(tolerance)
  n <- nrow(peaks(object))

  par_in$object <- object

  object <- do.call(set_general_bounds, par_in)

  for ( par in names(par_out) ) {
    if ( par %in% c("position", "width", "height", "fraction.gauss")) {
      expect_equal(lower_bounds(object)[[par]], 
                   rep(par_out[[par]][[1]], n))
      expect_equal(upper_bounds(object)[[par]], 
                   rep(par_out[[par]][[2]], n))
    }
  }

  object

}

#------------------------------------------------------------------------------
check_offset_bounds <- function(object, par_in, tolerance) {

  testthat_tolerance(tolerance)
  n <- nrow(peaks(object))

  par_in$object <- object

  object <- do.call(set_offset_bounds, par_in)

  for ( par in names(par_in) ) {
    if ( par %in% c("position", "width", "height", "fraction.gauss")) {

      lower <- par_in[[par]][[1]]
      upper <- par_in[[par]][[2]]

      if ( par_in$relative ) {
        lower <- lower*peaks(object)[[par]]
        upper <- upper*peaks(object)[[par]]
      } 

      expect_equal(lower_bounds(object)[[par]], 
                   peaks(object)[[par]] + lower)
      expect_equal(upper_bounds(object)[[par]], 
                   peaks(object)[[par]] + upper)
    }
  }

  object

}

#==============================================================================>
# Functions for generating nmr object

gen_resonance <- function() {

  id <- "R1"
  nmrresonance_1d("0.5 d 20", id = id)

}

gen_species <- function() {

  r1 <- nmrresonance_1d("0.2 d 20", id = "R1")
  r2 <- nmrresonance_1d("0.8 d 20", id = "R2")

  id <- "S1"
  nmrspecies_1d(list(r1, r2), id = id)

}

gen_data <- function() {

  object <- gen_resonance()

  n <- 200
  x <- seq(0, 1, length.out = n)
  y <- values(object, x)
  p <- tibble(direct.shift = x,
              intensity = y)

  d <- new("NMRData1D")
  d@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  d
}

gen_fit <- function() {

  object <- gen_species()
  nmrdata <- gen_data()

  nmrfit_1d(object, nmrdata, delay = TRUE)

}



################################################################################
# Resonance



#==============================================================================>
test_that("1d absolute general bounds work with resonance", {

  object <- gen_resonance()

  par_abs <- list(position = c(2, 3), width = c(0.005, 5), 
                  height = c(0, 10), fraction.gauss = c(0.1, 0.9),
                  widen = TRUE)

  object <- check_general_bounds(object, par_abs, par_abs, 1e-6)

})

#------------------------------------------------------------------------------
test_that("1d widen general bounds work with resonance", {

  nmroptions$sf <- 500

  object <- gen_resonance()

  par_abs <- list(position = c(0, 0.5), width = c(0, 0.5), 
                  height = c(0, 0.5), fraction.gauss = c(0.1, 0.9),
                  widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par_abs)

  par_wide <- list(position = c(-1, 1), width = c(-1, 1), 
                   height = c(-1, 1), fraction.gauss = c(0, 1),
                   widen = FALSE)

  object <- check_general_bounds(object, par_wide, par_abs, 1e-6)

})

#------------------------------------------------------------------------------
test_that("1d relative general bounds work with resonance", {

  nmroptions$sf <- 500

  object <- gen_resonance()

  p <- tibble(direct.shift = c(-1, 1), intensity = complex(re=c(-1, 1)))
  d <- new("NMRData1D")
  d@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  par_rel <- list(position = c(0, 1), width = c(0, 1), 
                  height = c(0, 1), nmrdata = d, widen = TRUE)

  par_out <- list(position = c(-1, 1), width = c(0, 2*nmroptions$sf), 
                  height = c(0, 2), nmrdata = d, widen = TRUE)

  object <- check_general_bounds(object, par_rel, par_out, 1e-6)

})

#==============================================================================>
test_that("1d absolute offset bounds work with resonance", {

  nmroptions$sf <- 500

  object <- gen_resonance()

  par <- list(position = c(-0.05, 0.05), width = c(-0.5, 1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = FALSE)

  object <- check_offset_bounds(object, par, 1e-6)

})

#------------------------------------------------------------------------------
test_that("1d relative offset bounds work with resonance", {

  nmroptions$sf <- 500

  object <- gen_resonance()

  par <- list(position = c(-0.01, 0.01), width = c(-0.1, 0.1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = TRUE)

  object <- check_offset_bounds(object, par, 1e-6)

})



################################################################################
# Species



#==============================================================================>
test_that("1d absolute general bounds work with species", {

  object <- gen_species()

  par_abs <- list(position = c(2, 3), width = c(0.005, 5), 
                  height = c(0, 10), fraction.gauss = c(0.1, 0.9),
                  widen = TRUE)

  object <- check_general_bounds(object, par_abs, par_abs, 1e-6)

})

#------------------------------------------------------------------------------
test_that("1d widen general bounds work with species", {

  object <- gen_species()

  par_abs <- list(position = c(0, 0.5), width = c(0, 0.5), 
                  height = c(0, 0.5), fraction.gauss = c(0.1, 0.9),
                  widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par_abs)

  par_wide <- list(position = c(-1, 1), width = c(-1, 1), 
                   height = c(-1, 1), fraction.gauss = c(0, 1),
                   widen = FALSE)

  object <- check_general_bounds(object, par_wide, par_abs, 1e-6)

})

#------------------------------------------------------------------------------
test_that("1d relative general bounds work with species", {

  object <- gen_species()

  p <- tibble(direct.shift = c(-1, 1), intensity = complex(re=c(-1, 1)))
  d <- new("NMRData1D")
  d@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  par_rel <- list(position = c(0, 1), width = c(0, 1), 
                  height = c(0, 1), nmrdata = d, widen = TRUE)

  par_out <- list(position = c(-1, 1), width = c(0, 2*nmroptions$sf), 
                  height = c(0, 2), nmrdata = d, widen = TRUE)

  object <- check_general_bounds(object, par_rel, par_out, 1e-6)

})

#==============================================================================>
test_that("1d absolute offset bounds work with species", {

  object <- gen_species()

  par <- list(position = c(-0.05, 0.05), width = c(-0.5, 1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = FALSE)

  object <- check_offset_bounds(object, par, 1e-6)

})

#------------------------------------------------------------------------------
test_that("1d relative offset bounds work with species", {

  object <- gen_species()

  par <- list(position = c(-0.01, 0.01), width = c(-0.1, 0.1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = TRUE)

  object <- check_offset_bounds(object, par, 1e-6)

})



################################################################################
# Fit



#==============================================================================>
test_that("1d absolute general bounds work with fit", {

  object <- gen_fit()

  par_abs <- list(position = c(2, 3), width = c(0.005, 5), 
                  height = c(0, 10), fraction.gauss = c(0.1, 0.9),
                  widen = TRUE)

  object <- check_general_bounds(object, par_abs, par_abs, 1e-6)

})

#------------------------------------------------------------------------------
test_that("1d widen general bounds work with fit", {

  object <- gen_fit()

  par_abs <- list(position = c(0, 0.5), width = c(0, 0.5), 
                  height = c(0, 0.5), fraction.gauss = c(0.1, 0.9),
                  widen = TRUE, object = object)

  object <- do.call(set_general_bounds, par_abs)

  par_wide <- list(position = c(-1, 1), width = c(-1, 1), 
                   height = c(-1, 1), fraction.gauss = c(0, 1),
                   widen = FALSE)

  object <- check_general_bounds(object, par_wide, par_abs, 1e-6)

})

#------------------------------------------------------------------------------
test_that("1d relative general bounds work with fit", {

  object <- gen_fit()

  p <- tibble(direct.shift = c(-1, 1), intensity = complex(re=c(-1, 1)))
  d <- new("NMRData1D")
  d@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  par_rel <- list(position = c(0, 1), width = c(0, 1), 
                  height = c(0, 1), nmrdata = d, widen = TRUE)

  par_out <- list(position = c(-1, 1), width = c(0, 2*nmroptions$sf), 
                  height = c(0, 2), nmrdata = d, widen = TRUE)

  object <- check_general_bounds(object, par_rel, par_out, 1e-6)

})

#==============================================================================>
test_that("1d absolute offset bounds work with fit", {

  object <- gen_fit()

  par <- list(position = c(-0.05, 0.05), width = c(-0.5, 1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = FALSE)

  object <- check_offset_bounds(object, par, 1e-6)

})

#------------------------------------------------------------------------------
test_that("1d relative offset bounds work with fit", {

  object <- gen_fit()

  par <- list(position = c(-0.01, 0.01), width = c(-0.1, 0.1), 
              height = c(-0.5, 0.5), widen = TRUE, relative = TRUE)

  object <- check_offset_bounds(object, par, 1e-6)

})
