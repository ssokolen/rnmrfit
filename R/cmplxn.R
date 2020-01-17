# Definition of a new set of complex classes that expand beyond the i notation


#' @importFrom vctrs vec_arith
NULL



#==============================================================================>
# 1D -- basically a copy of complex() to match the 2D version below
#==============================================================================>


#' @export
new_cmplx1 <- function(r = double(), i = double()) {
  vec_assert(r, ptype = double())
  vec_assert(i, ptype = double())
  
  new_rcrd(list(r = r, i = i), class = "vctrs_cmplx1")
}

#' @export
cmplx1 <- function(r = 0, i = 0) {
  c(r, i) %<-% vec_cast_common(r, i, .to = double())
  c(r, i) %<-% vec_recycle_common(r, i)
  
  new_cmplx1(r, i)
}

#' @method format vctrs_cmplx1
#' @export
format.vctrs_cmplx1 <- function(x, ...) {
  r <- field(x, "r")
  i <- field(x, "i")
  
  signs <- c("1"="+", "0"="+", "-1"="")
  f <- function (x) signs[as.character(sign(x))]
  out <- paste0(r, f(i), i, "i")
  out[is.na(r) | is.na(i)] <- NA
  
  out
}

#' @method vec_ptype_abbr vctrs_cmplx1
#' @export
vec_ptype_abbr.vctrs_cmplx1 <- function(x) "cmplx1"

#' @method vec_ptype_full vctrs_cmplx1
#' @export
vec_ptype_full.vctrs_cmplx1 <- function(x) "complex1d"

#' @method dim vctrs_cmplx1
#' @export
dim.vctrs_cmplx1 <- function(x) NULL

#------------------------------------------------------------------------------
# New subsetting functions

#' @export
`$.vctrs_cmplx1` <- function(x, name) field(x, name)

#' @export
`$<-.vctrs_cmplx1` <- function(x, name, value) {
  field(x, name) <- value
  x
}

#------------------------------------------------------------------------------
# Casting

#---------------------------------------
# Type definitions

#' @method vec_ptype2 vctrs_cmplx1
#' @export
#' @export vec_ptype2.vctrs_cmplx1
vec_ptype2.vctrs_cmplx1 <- function(x, y, ...) {
  UseMethod("vec_ptype2.vctrs_cmplx1", y)
}

#' @method vec_ptype2.vctrs_cmplx1 default
#' @export
vec_ptype2.vctrs_cmplx1.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  stop_incompatible_type(x, y, x_arg = x_arg, y_arg = y_arg)
}

#' @method vec_ptype2.vctrs_cmplx1 vctrs_unspecified
#' @export
vec_ptype2.vctrs_cmplx1.vctrs_unspecified <- function(x, y, ...) x

#' @method vec_ptype2.vctrs_cmplx1 vctrs_cmplx1
#' @export
vec_ptype2.vctrs_cmplx1.vctrs_cmplx1 <- function(x, y, ...) new_cmplx1()

#' @method vec_ptype2.vctrs_cmplx1 complex
#' @export
vec_ptype2.vctrs_cmplx1.complex <- function(x, y, ...) new_cmplx1()

#' @method vec_ptype2 complex
#' @export
#' @export vec_ptype2.complex
vec_ptype2.complex <- function(x, y, ...) {
  UseMethod("vec_ptype2.complex", y)
}

#' @method vec_ptype2.complex vctrs_cmplx1
#' @export
vec_ptype2.complex.vctrs_cmplx1 <- function(x, y, ...) new_cmplx1()

#' @method vec_ptype2.vctrs_cmplx1 double
#' @export
vec_ptype2.vctrs_cmplx1.double <- function(x, y, ...) new_cmplx1()

#' @method vec_ptype2.double vctrs_cmplx1
#' @export
vec_ptype2.double.vctrs_cmplx1 <- function(x, y, ...) new_cmplx1()

#---------------------------------------
# Cast definitions
vec_cast.vctrs_cmplx1 <- function(x, to) UseMethod("vec_cast.vctrs_cmplx1")
vec_cast.vctrs_cmplx1.default <- function(x, to) vec_default_cast(x, to)

vec_cast.vctrs_cmplx1.vctrs_cmplx1 <- function(x, to) x
vec_cast.numeric.vctrs_cmplx1 <- function(x, to) Re(x)
vec_cast.complex.vctrs_cmplx1 <- function(x, to) {
  complex(real = Re(x), imag = Im(x))
}
vec_cast.vctrs_cmplx1.numeric <- function(x, to) cmplx1(r = x)
vec_cast.vctrs_cmplx1.complex <- function(x, to) cmplx1(r = Re(x), i = Im(x))

#------------------------------------------------------------------------------
# Traditional Re()/Im()/Conj()

#' @export
Re.vctrs_cmplx1 <- function(z) z$r

#' @export
Im.vctrs_cmplx1 <- function(z) z$i

#' @export
Conj.vctrs_cmplx1 <- function(z) new_cmplx1(r = z$r, i = -z$i)

#------------------------------------------------------------------------------
# Summary

#' @export
as_tibble.vctrs_cmplx1 <- function(x, ...) as_tibble(unclass(x))

#' @export
summary.vctrs_cmplx1 <- function(object, ..., 
                                 digits = max(3, getOption("digits") - 3)) {
  summary(as_tibble(object))

}

#------------------------------------------------------------------------------
# Arithmetic

#---------------------------------------
# Boilerplate

#' @export
vec_arith.vctrs_cmplx1 <- function(op, x, y, ...) {
  UseMethod("vec_arith.vctrs_cmplx1", y)
}

#' @export
vec_arith.vctrs_cmplx1.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
vec_arith.complex <- function(op, x, y, ...) {
  UseMethod("vec_arith.complex", y)
}

#' @export
vec_arith.complex.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#---------------------------------------
# cmplx1

#' @S3method vec_arith.vctrs_cmplx1 vctrs_cmplx1
vec_arith.vctrs_cmplx1.vctrs_cmplx1 <- function(op, x, y, ...) {
  switch(
    op,
    "+" = new_cmplx1(r = x$r + y$r, i = x$i + y$i), 
    "-" = new_cmplx1(r = x$r - y$r, i = x$i - y$i), 
    "*" = new_cmplx1(r = x$r*y$r - x$i*y$i, i = x$r*y$i + x$i*y$r),
    "/" = (x *  Conj(y)) / (y * Conj(y))$r,
    stop_incompatible_op(op, x, y)
  )
}

#---------------------------------------
# numeric

#' @S3method vec_arith.vctrs_cmplx1 numeric
vec_arith.vctrs_cmplx1.numeric <- function(op, x, y, ...) {
  switch(
    op,
    "/" = new_cmplx1(r = x$r/y, i = x$i/y),
    vec_arith(op, x, cmplx1(r = y))
  )
}

#' @S3method vec_arith.numeric vctrs_cmplx1
vec_arith.numeric.vctrs_cmplx1 <- function(op, x, y, ...) {
  vec_arith(op, cmplx1(r = x), y)
}

#---------------------------------------
# complex 

#' @S3method vec_arith.vctrs_cmplx1 complex
vec_arith.vctrs_cmplx1.complex <- function(op, x, y, ...) {
   vec_arith(op, x, vec_cast(y, new_cmplx1())) 
}

#' @S3method vec_arith.complex vctrs_cmplx1
vec_arith.complex.vctrs_cmplx1 <- function(op, x, y, ...) {
  vec_arith(op, vec_cast(x, complex()), y)
}



#==============================================================================>
# 2D 
#==============================================================================>


#' @export
new_cmplx2 <- function(rr = double(), ri = double(), 
                       ir = double(), ii = double()) {
  vec_assert(rr, ptype = double())
  vec_assert(ri, ptype = double())
  vec_assert(ir, ptype = double())
  vec_assert(ii, ptype = double())
  
  new_rcrd(list(rr = rr, ri = ri, ir = ir, ii = ii), class = "vctrs_cmplx2")
}

#' @export
cmplx2 <- function(rr = 0, ri = 0, ir = 0, ii = 0) {
  c(rr, ri, ir, ii) %<-% vec_cast_common(rr, ri, ir, ii, .to = double())
  c(rr, ri, ir, ii) %<-% vec_recycle_common(rr, ri, ir, ii)
  
  new_cmplx2(rr, ri, ir, ii)
}

#' @export
format.vctrs_cmplx2 <- function(x, ...) {
  rr <- field(x, "rr")
  ri <- field(x, "ri")
  ir <- field(x, "ir")
  ii <- field(x, "ii")
  
  signs <- c("1"="+", "0"="+", "-1"="")
  f <- function (x) signs[as.character(sign(x))]
  out <- paste0(rr, f(ri), ri, "j", f(ir), ir, "i", f(ii), ii, "ji")
  out[is.na(rr) | is.na(ri) | is.na(ir) | is.na(ii)] <- NA
  
  out
}

#' @export
vec_ptype_abbr.vctrs_cmplx2 <- function(x) "cmplx2"

#' @export
vec_ptype_full.vctrs_cmplx2 <- function(x) "complex2d"

#' @export
vec_type2.vctrs_cmplx2.vctrs_unspecified <- function(x, y, ...) x

#------------------------------------------------------------------------------
# New subsetting functions

#' @export
`$.vctrs_cmplx2` <- function(x, name) field(x, name)
`$<-.vctrs_cmplx2` <- function(x, name, value) {
  field(x, name) <- value
  x
}

#------------------------------------------------------------------------------
# Casting

#---------------------------------------
# Type definitions
vec_type2.vctrs_cmplx2 <- function(x, y, ...) {
  UseMethod("vec_type2.vctrs_cmplx2", y)
}

vec_type2.vctrs_cmplx2.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  stop_incompatible_type(x, y, x_arg = x_arg, y_arg = y_arg)
}

vec_type2.vctrs_cmplx2.vctrs_unspecified <- function(x, y, ...) x
vec_type2.vctrs_cmplx2.vctrs_cmplx2 <- function(x, y, ...) new_cmplx2()

vec_type2.vctrs_cmplx2.complex <- function(x, y, ...) new_cmplx2()
vec_type2.complex.vctrs_cmplx2 <- function(x, y, ...) new_cmplx2()

vec_type2.vctrs_cmplx2.numeric <- function(x, y, ...) new_cmplx2()
vec_type2.numeric.vctrs_cmplx2 <- function(x, y, ...) new_cmplx2()

#---------------------------------------
# Cast definitions

vec_cast.vctrs_cmplx2 <- function(x, to) UseMethod("vec_cast.vctrs_cmplx2")
vec_cast.vctrs_cmplx2.default <- function(x, to) vec_default_cast(x, to)

vec_cast.vctrs_cmplx2.vctrs_cmplx2 <- function(x, to) x
vec_cast.numeric.vctrs_cmplx2 <- function(x, to) Re(x)
vec_cast.complex.vctrs_cmplx2 <- function(x, to) {
  complex(real = Re(x), imag = Im(x))
}
vec_cast.vctrs_cmplx2.numeric <- function(x, to) cmplx2(rr = x)
vec_cast.vctrs_cmplx2.complex <- function(x, to) cmplx2(rr = Re(x), ii = Im(x))

#------------------------------------------------------------------------------
# Traditional Re()/Im()

#' @export
Re.vctrs_cmplx2 <- function(z) z$rr

#' @export
Im.vctrs_cmplx2 <- function(z) z$ii

#------------------------------------------------------------------------------
# Summary

#' @export
as_tibble.vctrs_cmplx2 <- as_tibble.vctrs_cmplx1

#' @export
summary.vctrs_cmplx2 <- summary.vctrs_cmplx1
