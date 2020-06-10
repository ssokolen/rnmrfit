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

test_that("1d fit works from resonance", {

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
  d@processed <- p

  # Fitting
  r <- nmrresonance_1d('0.48 s')
  f <- nmrfit_1d(r, d)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_basic.htm")

  # Position is generally a good metric for convergence
  expect_true( abs(peaks(r_ideal)$position - peaks(f)$position) < 1e-3 )

  #----------------------------------------
  # Adding baseline
  p <- tibble(direct.shift = x,
              intensity = y + add_baseline(x, c(0, 0.2, 0.1)))

  d <- new("NMRData1D")
  d@processed <- p

  # Fitting
  r <- nmrresonance_1d('0.48 s')
  f <- nmrfit_1d(r, d)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_basic.htm")

  # Position is generally a good metric for convergence
  expect_true( abs(peaks(r_ideal)$position - peaks(f)$position) < 1e-3 )

  #----------------------------------------
  # Adding phase
  p <- tibble(direct.shift = x,
              intensity = y + add_baseline(x, c(0, 0.2, 0.1))) %>%
       mutate(intensity = add_phase(direct.shift, intensity, 0.1))

  d <- new("NMRData1D")
  d@processed <- p

  # Fitting
  r <- nmrresonance_1d('0.48 s')
  f <- nmrfit_1d(r, d)
  p <- plot(f)
  htmlwidgets::saveWidget(p, "1d_basic.htm")

  # Position is generally a good metric for convergence
  expect_true( abs(peaks(r_ideal)$position - peaks(f)$position) < 1e-3 )

  # Explicitly checking phase
  expect_true( abs(phase(f) + 0.1) < 0.01 )

})
