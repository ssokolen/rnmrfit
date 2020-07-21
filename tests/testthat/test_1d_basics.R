context("1D Basics")

#==============================================================================>
test_that("1d resonance construction works", {

  nmroptions$sf <- 500

  id <- "R1"
  resonance <- nmrresonance_1d("0.5 d 20", id = id)

  # Checking getters
  expect_equal(id(resonance), id)
  
  expect_equal(peaks(resonance), resonance@peaks)
  expect_true(all(peaks(resonance, TRUE)[,"resonance"] == id) ) 

  expect_equal(couplings(resonance), resonance@couplings)
  expect_true(all(couplings(resonance, TRUE)[,"resonance.1"] == id)) 
  expect_true(all(couplings(resonance, TRUE)[,"resonance.2"] == id)) 

})

#==============================================================================>
test_that("1d species construction works", {

  nmroptions$sf <- 500

  r1 <- nmrresonance_1d("0.2 d 20", id = "R1")
  r2 <- nmrresonance_1d("0.8 d 20", id = "R2")

  id <- "S1"
  species <- nmrspecies_1d(list(r1, r2), id = id)

  expect_equal(id(species), id)
  expect_equal(peaks(species), 
               rbind(peaks(r1, TRUE), peaks(r2, TRUE)))
  expect_equal(couplings(species), 
               rbind(couplings(r1, TRUE), couplings(r2, TRUE)))
  
  # Checking setters
  d <- peaks(r1)
  d$position = d$position + 0.1
  peaks(r1) <- d
  peaks(species) <- rbind(peaks(r1, TRUE), peaks(r2, TRUE))
  expect_equal(peaks(species), rbind(peaks(r1, TRUE), peaks(r2, TRUE)))

  expect_error(peaks(species) <- rbind(peaks(r1), peaks(r2)),
               'Input data.frame must have a "resonance" column.')

  r3 <- nmrresonance_1d("0.4 d 20", id = "R3")
  expect_warning(peaks(species) <- rbind(peaks(r1, TRUE), peaks(r3, TRUE)),
                 "The following resonance is not defined, ignoring: R3")

  # Checking conformity checks
  p <- tibble(direct.shift = c(-1, 1), intensity = complex(re=c(-1, 1)))
  d <- new("NMRData1D")
  d@acqus <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p
  expect_true(check_conformity(species, d))

  p$direct.shift <- c(-1, 0)
  d@processed <- p
  expect_error(check_conformity(species, d),
               "Some peaks fall outside the data's chemical shift")
})

#==============================================================================>
test_that("1d fit construction works", {

  nmroptions$sf <- 500

  #----------------------------------------
  # Basic peak
  r_ideal <- nmrresonance_1d('0.5 s', width = 3)
  d <- nmrdata_1d_from_scaffold(r_ideal)

  # Fitting
  r <- nmrresonance_1d('0.485 s')
  fit <- nmrfit_1d(r, d, delay = TRUE)

  # Checking getters
  expect_equal(peaks(r)$position, peaks(fit)$position)
})

