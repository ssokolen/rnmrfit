context("1D Fit")

# Add complex noise
add_noise <- function(y, sd) {
  n <- length(y)
  y + cmplx1(r = rnorm(n, 0, sd), i = rnorm(n, 0, sd))
}

add_baseline <- function(x, order) {
  basis <- splines::bs(x, degree = length(order) - 1, intercept = TRUE)
  baseline <- (basis %*% order)[, 1]
  cmplx1(r = baseline, i = baseline)
}

add_phase <- function(x, y, theta) {
  if ( length(theta) > 1 ) {
    theta = theta[1] + theta[2]*x
  }
  cmplx1(r = Re(y)*cos(theta) + Im(y)*sin(theta), 
         i = -Re(y)*sin(theta) + Im(y)*cos(theta))
}

#==============================================================================>
test_that("1d fit construction works", {

  nmroptions$sf <- 500

  #----------------------------------------
  # Basic peak
  r_ideal <- nmrresonance_1d('0.5 s', width = 3)

  n <- 200
  x <- seq(0, 1, length.out = n)
  y <- values(r_ideal, x) %>% add_noise(0.02)
  p <- tibble(direct.shift = x,
              intensity = y)

  d <- new("NMRData1D")
  d@procs <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  # Fitting
  r <- nmrresonance_1d('0.485 s')
  fit <- nmrfit_1d(r, d, delay = TRUE)

  # Checking absolute general bounds
  fit <- set_general_bounds(fit, 
    position = c(0, 5), width = c(0, 5), height = c(0, 5), 
  )
  expect_true(all(lower_bounds(fit)$position == 0))
  expect_true(all(lower_bounds(fit)$width == 0))
  expect_true(all(lower_bounds(fit)$heigft == 0))

  expect_true(all(upper_bounds(fit)$position == 5))
  expect_true(all(upper_bounds(fit)$width == 5))
  expect_true(all(upper_bounds(fit)$heigft == 5))

})


#==============================================================================>
test_that("1d fit works with singlets", {

  nmroptions$sf <- 500

  #----------------------------------------
  # Basic peak
  r_ideal <- nmrresonance_1d('0.5 s', width = 3)

  n <- 200
  x <- seq(0, 1, length.out = n)
  y <- values(r_ideal, x) %>% add_noise(0.02)
  p <- tibble(direct.shift = x,
              intensity = y)

  d <- new("NMRData1D")
  d@procs <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  # Fitting
  r <- nmrresonance_1d('0.485 s')
  f <- nmrfit_1d(r, d)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_singlet.htm")

  # Position is generally a good metric for convergence
  expect_true( abs(peaks(r_ideal)$position - peaks(f)$position) < 1e-3 )

  #----------------------------------------
  # Adding baseline
  p <- tibble(direct.shift = x,
              intensity = y + add_baseline(x, c(0, 0.2, 0.1)))

  d <- new("NMRData1D")
  d@procs <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  # Fitting
  f <- nmrfit_1d(r, d)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_singlet.htm")

  # Position is generally a good metric for convergence
  expect_true( abs(peaks(r_ideal)$position - peaks(f)$position) < 1e-3 )

  #----------------------------------------
  # Adding phase
  p <- tibble(direct.shift = x,
              intensity = y + add_baseline(x, c(0, 0.2, 0.1))) %>%
       mutate(intensity = add_phase(direct.shift, intensity, 0.1))

  d <- new("NMRData1D")
  d@procs <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  # Fitting
  f <- nmrfit_1d(r, d)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_singlet.htm")

  # Position is generally a good metric for convergence
  expect_true( abs(peaks(r_ideal)$position - peaks(f)$position) < 1e-3 )

  # Explicitly checking phase
  expect_true( abs(phase(f) + 0.1) < 0.03 )

})

#==============================================================================>
test_that("1d fit works with multiplets", {

  nmroptions$sf <- 500

  #----------------------------------------
  # Basic peak
  r_ideal <- nmrresonance_1d('0.5 t 15', width = 3)

  n <- 200
  x <- seq(0, 1, length.out = n)
  y <- values(r_ideal, x) %>% add_noise(0.02)
  p <- tibble(direct.shift = x,
              intensity = y)

  d <- new("NMRData1D")
  d@procs <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  # Fitting
  r <- nmrresonance_1d('0.485 t 20', position.leeway = 0.5)
  f <- nmrfit_1d(r, d)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_multiplet.htm")

  # Position is generally a good metric for convergence
  expect_true( all(abs(peaks(r_ideal)$position - peaks(f)$position) < 2e-3) )

})

#==============================================================================>
test_that("1d fit works with species", {

  nmroptions$sf <- 500

  #----------------------------------------
  # Basic peak
  r1 <- nmrresonance_1d('0.2 t 15', width = 3)
  r2 <- nmrresonance_1d('0.8 d 25', width = 3)
  peaks(r2) <- peaks(r2) %>% mutate(height = height*2)
  s_ideal <- nmrspecies_1d(list(r1, r2))

  n <- 200
  x <- seq(0, 1, length.out = n)
  y <- values(s_ideal, x) %>% add_noise(0.02)
  p <- tibble(direct.shift = x,
              intensity = y)

  d <- new("NMRData1D")
  d@procs <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  # Fitting
  r1 <- nmrresonance_1d('0.185 t 20', position.leeway = 0.5)
  r2 <- nmrresonance_1d('0.785 d 20', position.leeway = 0.5)
  s <- nmrspecies_1d(list(r1, r2), areas = c(1, 3), connections.leeway = 0.5)

  f <- nmrfit_1d(s, d)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_species.htm")

  # Position is generally a good metric for convergence
  expect_true( all(abs(peaks(s_ideal)$position - peaks(f)$position) < 1e-3) )

  #----------------------------------------
  # Adding baseline
  p <- tibble(direct.shift = x,
              intensity = y + add_baseline(x, c(0, 0.2, 0.1)))

  d <- new("NMRData1D")
  d@procs <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  # Fitting
  f <- nmrfit_1d(s, d)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_species.htm")

  # Position is generally a good metric for convergence
  expect_true( all(abs(peaks(s_ideal)$position - peaks(f)$position) < 1e-3) )

  #----------------------------------------
  # Adding phase
  p <- tibble(direct.shift = x,
              intensity = y + add_baseline(x, c(0, 0.2, 0.1))) %>%
       mutate(intensity = add_phase(direct.shift, intensity, c(0.1, 0.2)))

  d <- new("NMRData1D")
  d@procs <- list(direct = list(sfo1 = nmroptions$sf))
  d@processed <- p

  # Fitting
  f <- nmrfit_1d(s, d, phase.order = 1)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_species.htm")

  # Position is generally a good metric for convergence
  expect_true( all(abs(peaks(s_ideal)$position - peaks(f)$position) < 1e-3) )

  # Explicitly checking phase
  expect_true( abs(phase(f)[1] + 0.1) < 0.02 )

})

