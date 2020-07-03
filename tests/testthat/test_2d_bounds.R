context("2D Bounds")

nmroptions$sf <- 500

#==============================================================================>
# The same tests wll be run on resonance/species/fit
check_general_bounds <- function(object, par, tolerance) {

  testthat_tolerance(tolerance)
  n <- nrow(peaks(object))

  for ( p in names(par) ) {
    if ( p %in% c("position", "width", "height", "fraction.gauss")) {

      expect_equal(lower_bounds(object)[[p]], rep(par[[p]][[1]], n))
      expect_equal(upper_bounds(object)[[p]], rep(par[[p]][[2]], n))
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
  r1.direct <- nmrresonance_1d("0.5 d 20", id = id)
  r1.indirect <- nmrresonance_1d("0.7 d 20", id = id)
  nmrresonance_2d(r1.direct, r1.indirect)

}

gen_species <- function() {

  id <- "R1"
  r1.direct <- nmrresonance_1d("0.5 d 20", id = id)
  r1.indirect <- nmrresonance_1d("0.7 d 20", id = id)
  r1 <- nmrresonance_2d(r1.direct, r1.indirect)

  id <- "R2"
  r2.direct <- nmrresonance_1d("0.2 s", id = id)
  r2.indirect <- nmrresonance_1d("0.3 d 20", id = id)
  r2 <- nmrresonance_2d(r2.direct, r2.indirect)

  id <- "S1"
  nmrspecies_2d(list(r1, r2), id = id)

}

gen_data <- function() {

  x1 <- seq(0, 1, length.out = 50)
  x2 <- seq(0, 1, length.out = 50)
  p <- expand.grid(x1 = x1, x2 = x2)

  intensity <- values(species, p$x1, p$x2, components = 'rr/ri/ir/ii')

  p <- tibble(direct.shift = p$x1, indirect.shift = p$x2,
              intensity = intensity)

  d <- new("NMRData2D")
  d@processed <- p

  d
}

gen_fit <- function() {

  object <- gen_species()
  nmrdata <- gen_data()

  nmrfit_2d(object, nmrdata, delay = TRUE)

}



################################################################################
# Resonance



#==============================================================================>
test_that("1d absolute general bounds work with resonance", {

  par.abs <- list(position = c(2, 3), width = c(0.005, 5), 
                  height = c(0, 10), fraction.gauss = c(0.1, 0.9),
                  widen = TRUE)

  object <- gen_resonance()
  opts <- par.abs
  opts$object <- object
  object <- do.call(set_general_bounds, opts)

  check_general_bounds(direct(object), par.abs, 1e-6)
  check_general_bounds(indirect(object), par.abs, 1e-6)

})

if ( FALSE ) {
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
}
