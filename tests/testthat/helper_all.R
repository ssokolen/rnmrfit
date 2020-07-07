#==============================================================================>
# Preset 1D objects

#---------------------------------------
gen_resonance_1d <- function() {

  id <- "R1"
  nmrresonance_1d("0.5 d 20", id = id)
}

#---------------------------------------
gen_species_1d <- function() {

  r1 <- nmrresonance_1d("0.3 d 20", id = "R1")
  r2 <- nmrresonance_1d("0.7 d 20", id = "R2")

  id <- "S1"
  nmrspecies_1d(list(r1, r2), id = id)
}

#---------------------------------------
gen_data_1d <- function(object) {

  n <- 200

  x.min <- min(peaks(object)$position) - 0.2
  x.max <- max(peaks(object)$position) + 0.2

  x <- seq(x.min, x.max, length.out = n)
  y <- values(object, x)
  p <- tibble(direct.shift = x, intensity = y)

  d <- new("NMRData1D")
  d@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  d
}

#---------------------------------------
gen_fit_1d <- function(object) {

  nmrdata <- gen_data_1d(object)

  nmrfit_1d(object, nmrdata, delay = TRUE)
}

#==============================================================================>
# Preset 2D objects

#---------------------------------------
gen_resonance_2d <- function() {

  id <- "R1"
  r1.direct <- nmrresonance_1d("0.5 d 20", id = id)
  r1.indirect <- nmrresonance_1d("0.7 d 20", id = id)
  nmrresonance_2d(r1.direct, r1.indirect)
}

#---------------------------------------
gen_species_2d <- function() {

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

#---------------------------------------
gen_data_2d <- function(object) {

  n <- 50

  x.min <- min(peaks(direct(object))$position) - 0.2
  x.max <- max(peaks(direct(object))$position) + 0.2
  x1 <- seq(x.min, x.max, length.out = n)

  x.min <- min(peaks(indirect(object))$position) - 0.2
  x.max <- max(peaks(indirect(object))$position) + 0.2
  x2 <- seq(x.min, x.max, length.out = n)
  
  p <- expand.grid(x1 = x1, x2 = x2)
  intensity <- values(object, p$x1, p$x2, components = 'rr/ri/ir/ii')

  p <- tibble(direct.shift = p$x1, indirect.shift = p$x2,
              intensity = intensity)

  d <- new("NMRData2D")
  d@processed <- p

  d
}

#---------------------------------------
gen_fit_2d <- function(object) {

  nmrdata <- gen_data_2d(object)

  nmrfit_2d(object, nmrdata, delay = TRUE)
}

