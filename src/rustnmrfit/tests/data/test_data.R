library(jsonlite)
library(numDeriv)
library(RcppFaddeeva)
library(splines)

#==============================================================================>
# Lineshape functions

#---------------------------------------
f_lorentz <- function(x, p) {
  z = (p[1] - x)/p[2] 
  p[3]*complex(re = 1, im = z)/(z*z + 1)
}

#---------------------------------------
f_voigt <- function(x, p) {
  wg <- p[2] * p[4] / (1 - p[4])

  zn <- complex(re = 0, im = p[2]/(sqrt(2)*wg))
  fn <- Faddeeva_w(zn)

  z <- complex(re = p[1] - x, im = p[2])/complex(re = (sqrt(2)*wg), 0)
  p[3]*Faddeeva_w(z)/fn
}

#---------------------------------------
f_normalize <- function(y) {
  (y - min(y))/(max(y) - min(y))*2-1
}

#---------------------------------------
f_phase <- function(x, y, theta) {
  theta <- theta[1] + theta[2]*x
  complex(re = Re(y)*cos(theta) + Im(y)*sin(theta),
          im = -Re(y)*sin(theta) + Im(y)*cos(theta))
}

#==============================================================================>
# Testing lineshapes

gen_data <- function(x, p, f, filename) {

  out <- list(par = p, x = x)

  # Generating y data
  y <- f(x, p)
  out$y_re <- Re(y)
  out$y_im <- Im(y)

  dy_re <- sapply(x, 
    function(x) {
      f_inner <- function(p) { Re(f(x, p)) }
      grad(f_inner, p) 
    })
  out$dy_re <- as.vector(t(dy_re))

  dy_im <- sapply(x, 
    function(x) {
      f_inner <- function(p) { Im(f(x, p)) }
      grad(f_inner, p) 
    })
  out$dy_im <- as.vector(t(dy_im))

  json <- toJSON(out, digits = 12)

  con <- file(filename)
  writeLines(json, con)
  close(con)
}

#------------------------------------------------------------------------------

x <- seq(0, 1, length.out = 100)

# Lorentz
p1 <- c(0.6, 0.1, 0.8, 0.0)
gen_data(x, p1, f_lorentz, "lorentz_data.json")

# Voigt
p2 <- c(0.3, 0.1, 0.6, 0.5)
gen_data(x, p2, f_voigt, "voigt_data.json")

# Combined
p3 <- c(p1, p2)
f <- function(x, p) f_lorentz(x, p[1:4]) + f_voigt(x, p[5:8])
gen_data(x, p3, f, "combined_data.json")

#==============================================================================>
# Testing fit

gen_data <- function(x, p, f, nb = 0, np = 0, filename) {

  out <- list(par = c(p, rep(0, nb*2 + np)), x = x)

  # Generating y data
  y <- f(x, p)
  out$y_re <- Re(y)
  out$y_im <- Im(y)

  # Generating basis
  if ( nb == 0 ) {
    basis <- matrix(0, nrow = length(x), ncol = 1)
  } else {
    basis <- bs(x, degree = nb - 1, intercept = TRUE)
  }
  out$basis <- as.vector(t(basis))

  # Adding length terms
  out$nb <- nb
  out$np <- np

  json <- toJSON(out, digits = 12)

  con <- file(filename)
  writeLines(json, con)
  close(con)
}

#------------------------------------------------------------------------------

set.seed(1111)
x <- seq(0, 1, length.out = 100)

# Lorentz
p1 <- c(0.6, 0.1, 0.8, 0.0)
f <- function(x, p) f_lorentz(x, p) + rnorm(length(x), 0, p1[3]/50)
gen_data(x, p1, f, 0, 0, "lorentz_fit.json")

# Voigt
p2 <- c(0.3, 0.05, 0.6, 0.5)
f <- function(x, p) f_voigt(x, p) + rnorm(length(x), 0, p2[3]/50)
gen_data(x, p2, f, 0, 0, "voigt_fit.json")

# Combined
p3 <- c(p1, p2)
f <- function(x, p) {
  f_lorentz(x, p[1:4]) + f_voigt(x, p[5:8]) + rnorm(length(x), 0, p3[3]/50)
}
gen_data(x, p3, f, 0, 0, "combined_fit.json")

# Baseline
p3 <- c(p1, p2)
f <- function(x, p) {
  f_lorentz(x, p[1:4]) + f_voigt(x, p[5:8]) + 
    rnorm(length(x), 0, p3[3]/50) + 
    f_normalize(-(x-.4)*(x-.4))*0.2
}
gen_data(x, p3, f, 3, 0, "baseline_fit.json")

# Phase
p3 <- c(p1, p2)
f <- function(x, p) {
  y <- f_lorentz(x, p[1:4]) + f_voigt(x, p[5:8]) +
    rnorm(length(x), 0, p3[3]/50) + 
    f_normalize(-(x-.4)*(x-.4))*0.2
  y <- f_phase(x, y, c(pi/12, pi/16))
}

gen_data(x, p3, f, 3, 2, "phase_fit.json")
