context("2D Basics")

#==============================================================================>
test_that("2d resonance construction works", {

  nmroptions$sf <- 500

  id <- "R1"
  r.direct <- nmrresonance_1d("0.5 d 20", id = id)
  r.indirect <- nmrresonance_1d("0.7 d 20", id = id)

  resonance <- nmrresonance_2d(r.direct, r.indirect)

  # Checking getters
  expect_equal(id(resonance), id)
  
  expect_equal(peaks(resonance), 
               rbind(cbind(dimension = "direct", peaks(r.direct)),
                     cbind(dimension = "indirect", peaks(r.indirect))))
  expect_true(all(peaks(resonance, TRUE)[,"resonance"] == id) ) 

  expect_equal(couplings(resonance),
               rbind(cbind(dimension = "direct", couplings(r.direct)),
                     cbind(dimension = "indirect", couplings(r.indirect))))
  expect_true(all(couplings(resonance, TRUE)[,"resonance.1"] == id)) 
  expect_true(all(couplings(resonance, TRUE)[,"resonance.2"] == id)) 

})

#==============================================================================>
test_that("2d species construction works", {

  nmroptions$sf <- 500

  id <- "R1"
  r.direct <- nmrresonance_1d("0.5 d 20", id = id)
  r.indirect <- nmrresonance_1d("0.7 d 20", id = id)
  r1 <- nmrresonance_2d(r.direct, r.indirect)

  id <- "R2"
  r.direct <- nmrresonance_1d("0.2 s", id = id)
  r.indirect <- nmrresonance_1d("0.3 d 20", id = id)
  r2 <- nmrresonance_2d(r.direct, r.indirect)

  id <- "S1"
  species <- nmrspecies_2d(list(r1, r2), id = id)

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

  r3 <- nmrresonance_2d(r.direct, r.indirect, id = "R3")

  expect_warning(peaks(species) <- rbind(peaks(r1, TRUE), peaks(r3, TRUE)),
                 "The following resonance is not defined, ignoring: R3")

  expect_error(peaks(species) <- peaks(r.direct, TRUE),
               'Input data.frame must have a "dimension" column.')

})

#==============================================================================>
test_that("2d projection works", {

  nmroptions$sf <- 500

  id <- "R1"
  r1.direct <- nmrresonance_1d("0.5 d 20", id = id)
  r1.indirect <- nmrresonance_1d("0.7 d 20", id = id)
  r1 <- nmrresonance_2d(r1.direct, r1.indirect)

  id <- "R2"
  r2.direct <- nmrresonance_1d("0.2 s", id = id)
  r2.indirect <- nmrresonance_1d("0.3 d 20", id = id)
  r2 <- nmrresonance_2d(r2.direct, r2.indirect)

  id <- "S1"
  species <- nmrspecies_2d(list(r1, r2), id = id)
  species.direct <- nmrspecies_1d(list(r1.direct, r2.direct), id = id)
  species.indirect <- nmrspecies_1d(list(r1.indirect, r2.indirect), id = id)

  expect_equal(peaks(direct(species)), peaks(species.direct))
  expect_equal(peaks(indirect(species)), peaks(species.indirect))

})
