library(MASS)
library(numDeriv)
library(RcppFaddeeva)

#---------------------------------------
# General Voigt function
f_voigt <- function(x, p) {
  wg <- p[2] * p[4] / (1 - p[4])

  zn <- complex(re = 0, im = p[2]/(sqrt(2)*wg))
  fn <- Faddeeva_w(zn)

  z <- complex(re = p[1] - x, im = p[2])/complex(re = (sqrt(2)*wg), 0)
  p[3]*Faddeeva_w(z)/fn
}

# Area of Voigt function in the real domain
f_voigt_area <- function(p) {

  wg <- p[2] * p[4] / (1 - p[4])

  zn <- complex(re = 0, im = p[2]/(sqrt(2)*wg))
  fn <- Faddeeva_w(zn)

  Re(sqrt(2*pi) * wg * p[3] / fn)

}

# Area of Voigt function in the real domain
f_voigt_area_deriv <- function(p) {

  wg <- p[2] * p[4] / (1 - p[4])

  zn <- complex(re = 0, im = p[2]/(sqrt(2)*wg))
  fn <- Faddeeva_w(zn)

  a = sqrt(2*pi) /  Re(fn);
  b = a/pi
  w = p[2]
  h = p[3]
  f = p[4]

  c(0, 
    a * h * f/(1-f), 
    a * wg,
    a * w * h * ( 1/( (1-f)^2 ) + 1/( f^2 ) - b/( f*(1-f) ) )
  ) 
}

#---------------------------------------
# Check that area function is correct

p <- runif(4)
print(p)

cat("\nArea comparison\n\n")
cat(sprintf('Analytical: %.4f. Numerical: %.4f\n',
            f_voigt_area(p), 
            integrate(function(x) {Re(f_voigt(x, p))}, -Inf, Inf)$value))

#---------------------------------------
# Check area derivates

cat("\nArea derivative comparison\n\n")
cat(sprintf('Analytical: %.4f. Numerical: %.4f\n',
            f_voigt_area_deriv(p), 
            grad(function(x) {Re(f_voigt_area(x))}, p)))

