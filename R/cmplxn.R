# Definition of a new set of complex classes that expand beyond the i notation

#' @import vctrs
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

#' @method vec_ptype2 complex
#' @export
#' @export vec_ptype2.complex
vec_ptype2.complex <- function(x, y, ...) {
  UseMethod("vec_ptype2.complex", y)
}

#' @method vec_ptype2.vctrs_cmplx1 default
#' @export
vec_ptype2.vctrs_cmplx1.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  stop_incompatible_type(x, y, x_arg = x_arg, y_arg = y_arg)
}

#' @method vec_ptype2.complex default
#' @export
vec_ptype2.complex.default <- function(x, y, ..., x_arg = "", y_arg = "") {
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

#' @method vec_cast vctrs_cmplx1
#' @export
#' @export vec_cast.vctrs_cmplx1
vec_cast.vctrs_cmplx1 <- function(x, to, ...) UseMethod("vec_cast.vctrs_cmplx1")

#' @method vec_cast.vctrs_cmplx1 default
#' @export
vec_cast.vctrs_cmplx1.default <- function(x, to, ...) vec_default_cast(x, to)

#' @method vec_cast.vctrs_cmplx1 vctrs_cmplx1
#' @export
vec_cast.vctrs_cmplx1.vctrs_cmplx1 <- function(x, to, ...) x

#' @method vec_cast.vctrs_cmplx1 complex
#' @export
vec_cast.vctrs_cmplx1.complex <- function(x, to, ...) {
  cmplx1(r = Re(x), i = Im(x))
}

#' @method vec_cast.complex vctrs_cmplx1
#' @export
vec_cast.complex.vctrs_cmplx1 <- function(x, to, ...) {
  complex(real = Re(x), imag = Im(x))
}

#' @method vec_cast.vctrs_cmplx1 double 
#' @export
vec_cast.vctrs_cmplx1.double <- function(x, to, ...) cmplx1(r = x)

#' @method vec_cast.double vctrs_cmplx1
#' @export
vec_cast.double.vctrs_cmplx1 <- function(x, to, ...) Re(x)


#------------------------------------------------------------------------------
# Traditional Re()/Im()/Conj()

#' @method Re vctrs_cmplx1
#' @export
Re.vctrs_cmplx1 <- function(z) z$r

#' @method Im vctrs_cmplx1
#' @export
Im.vctrs_cmplx1 <- function(z) z$i

#' @method Conj vctrs_cmplx1
#' @export
Conj.vctrs_cmplx1 <- function(z) new_cmplx1(r = z$r, i = -z$i)

#------------------------------------------------------------------------------
#' Flexible domain selection
#' 
#' @param x cmplx1 or cmplx2 object.
#' @param domain A character specifying which real/imaginary domains to output
#'               for intensity data. For 1D data "r/i" outputs both real and
#'               imaginary data, with "r" and "i" used to output one of either
#'               real or imaginary. For 2D data, this is expanded to
#'               "rr/ri/ir/ii". All domains are outputted by default. Warning:
#'               this function can be used to flip real and imaginary data by
#'               providing i/r as opposed to r/i.
#' @param use.cmplx1 cmplx1 is the internal representation of complex() that
#'                   matches cmplx2 for 2D data. It is replaced by default to
#'                   complex() for user-facing functions. However, this behavior
#'                   can be overriden here.
#' 
#' @return If one domain selected, a simple numeric vector. If two selected, a
#'         complex() object with real and imaginary term replaced by the
#'         selected variants. If 3 or 4 selected, a cmplx2 object.
#' 
#' @name domain
#' @export
domain <- function (x, ...) {
  UseMethod("domain", x)
}

#' @rdname domain
#' @export
domain.vctrs_cmplx1 <- function(x, domain = "r/i", use.cmplx1 = FALSE) {

  domain <- str_split(domain, "[ \t/_-]+")[[1]]

  err <- '"domain" must be composed of "r" or "i" characters split by /'

  if (! all(domain %in% c("r", "i")) ) stop(err)

  if ( length(domain) == 1 ) {
    get(domain, x)
  } else if ( length(domain) == 2) {
    if ( use.cmplx1 ) {
      cmplx1(r = get(domain[1], x), i = get(domain[2], x)) 
    } else {
      complex(re = get(domain[1], x), im = get(domain[2], x))
    }
  } else {
    stop(err)
  }
}

#------------------------------------------------------------------------------
# Summary

#' @method as_tibble vctrs_cmplx1
#' @export
as_tibble.vctrs_cmplx1 <- function(x, ...) as_tibble(unclass(x))

#' @method summary vctrs_cmplx1
#' @export
summary.vctrs_cmplx1 <- function(object, ..., 
                                 digits = max(3, getOption("digits") - 3)) {
  summary(as_tibble(object))

}

#------------------------------------------------------------------------------
# Arithmetic

#---------------------------------------
# Boilerplate

#' @method vec_arith vctrs_cmplx1
#' @export
#' @export vec_arith.vctrs_cmplx1
vec_arith.vctrs_cmplx1 <- function(op, x, y, ...) {
  UseMethod("vec_arith.vctrs_cmplx1", y)
}

#' @method vec_arith.vctrs_cmplx1 default
#' @export
vec_arith.vctrs_cmplx1.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @method vec_arith complex
#' @export
#' @export vec_arith complex
vec_arith.complex <- function(op, x, y, ...) {
  UseMethod("vec_arith.complex", y)
}

#' @method vec_arith.complex default
#' @export
vec_arith.complex.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @method vec_arith double
#' @export
#' @export vec_arith double
vec_arith.double <- function(op, x, y, ...) {
  UseMethod("vec_arith.double", y)
}

#' @method vec_arith.double default
#' @export
vec_arith.double.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#---------------------------------------
# cmplx1

#' @method vec_arith.vctrs_cmplx1 vctrs_cmplx1
#' @export
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
# complex 

#' @method vec_arith.vctrs_cmplx1 complex
#' @export
vec_arith.vctrs_cmplx1.complex <- function(op, x, y, ...) {
   vec_arith(op, x, vec_cast(y, new_cmplx1())) 
}

#' @method vec_arith.complex vctrs_cmplx1
#' @export
vec_arith.complex.vctrs_cmplx1 <- function(op, x, y, ...) {
  vec_arith(op, vec_cast(x, new_cmplx1()), y)
}

#---------------------------------------
# double

#' @method vec_arith.vctrs_cmplx1 double
#' @export
vec_arith.vctrs_cmplx1.double <- function(op, x, y, ...) {
  switch(
    op,
    "/" = new_cmplx1(r = x$r/y, i = x$i/y),
    vec_arith(op, x, cmplx1(r = y))
  )
}

#' @method vec_arith.double vctrs_cmplx1
#' @export
vec_arith.double.vctrs_cmplx1 <- function(op, x, y, ...) {
  vec_arith(op, cmplx1(r = x), y)
}



#------------------------------------------------------------------------------
# Math

#' @method vec_math vctrs_cmplx1
#' @export
vec_math.vctrs_cmplx1 <- function(.fn, .x, ...) {
  switch(.fn,
    sum = cmplx1(r = sum(.x$r), i = sum(.x$i)),
    mean = cmplx1(r = sum(.x$r), i = sum(.x$i))/as.numeric(length(.x)),
    vec_math_base(.fn, .x, ...)
  )
}

methods::setOldClass(c("vctrs_cmplx1", "vctrs_vctr"))



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

#' @method vec_ptype2 vctrs_cmplx2
#' @export
#' @export vec_ptype2.vctrs_cmplx2
vec_ptype2.vctrs_cmplx2 <- function(x, y, ...) {
  UseMethod("vec_ptype2.vctrs_cmplx2", y)
}

#' @method vec_ptype2.vctrs_cmplx2 default
#' @export
vec_ptype2.vctrs_cmplx2.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  stop_incompatible_type(x, y, x_arg = x_arg, y_arg = y_arg)
}

#' @method vec_ptype2.complex default
#' @export
vec_ptype2.complex.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  stop_incompatible_type(x, y, x_arg = x_arg, y_arg = y_arg)
}

#' @method vec_ptype2.vctrs_cmplx2 vctrs_unspecified
#' @export
vec_ptype2.vctrs_cmplx2.vctrs_unspecified <- function(x, y, ...) x

#' @method vec_ptype2.vctrs_cmplx2 vctrs_cmplx2
#' @export
vec_ptype2.vctrs_cmplx2.vctrs_cmplx2 <- function(x, y, ...) new_cmplx2()

#' @method vec_ptype2.vctrs_cmplx2 complex
#' @export
vec_ptype2.vctrs_cmplx2.complex <- function(x, y, ...) new_cmplx2()

#' @method vec_ptype2.complex vctrs_cmplx2
#' @export
vec_ptype2.complex.vctrs_cmplx2 <- function(x, y, ...) new_cmplx2()

#' @method vec_ptype2.vctrs_cmplx2 double
#' @export
vec_ptype2.vctrs_cmplx2.double <- function(x, y, ...) new_cmplx2()

#' @method vec_ptype2.double vctrs_cmplx2
#' @export
vec_ptype2.double.vctrs_cmplx2 <- function(x, y, ...) new_cmplx2()

#---------------------------------------
# Cast definitions

#' @method vec_cast vctrs_cmplx2
#' @export
#' @export vec_cast.vctrs_cmplx2
vec_cast.vctrs_cmplx2 <- function(x, to, ...) UseMethod("vec_cast.vctrs_cmplx2")

#' @method vec_cast.vctrs_cmplx2 default
#' @export
vec_cast.vctrs_cmplx2.default <- function(x, to, ...) vec_default_cast(x, to)

#' @method vec_cast.vctrs_cmplx2 vctrs_cmplx2
#' @export
vec_cast.vctrs_cmplx2.vctrs_cmplx2 <- function(x, to, ...) x

#' @method vec_cast.vctrs_cmplx2 complex
#' @export
vec_cast.vctrs_cmplx2.complex <- function(x, to, ...) {
  cmplx2(rr = Re(x), ii = Im(x))
}

#' @method vec_cast.complex vctrs_cmplx2
#' @export
vec_cast.complex.vctrs_cmplx2 <- function(x, to, ...) {
  complex(real = Re(x), imag = Im(x))
}

#' @method vec_cast.vctrs_cmplx2 double 
#' @export
vec_cast.vctrs_cmplx2.double <- function(x, to, ...) cmplx2(rr = x)

#' @method vec_cast.double vctrs_cmplx2
#' @export
vec_cast.double.vctrs_cmplx2 <- function(x, to, ...) Re(x)


#------------------------------------------------------------------------------
# Traditional Re()/Im()/Conj()

#' @method Re vctrs_cmplx2
#' @export
Re.vctrs_cmplx2 <- function(z) z$rr

#' @method Im vctrs_cmplx2
#' @export
Im.vctrs_cmplx2 <- function(z) z$ii

#------------------------------------------------------------------------------
#' Flexible domain selection

#' @rdname domain
#' @export
domain.vctrs_cmplx2 <- function(x, domain = "rr/ri/ir/ii", use.cmplx1 = FALSE) {

  domain <- str_split(domain, "[ \t/_-]+")[[1]]

  err <- '"domain" must be composed of "rr","ri","ir","ii" characters split by /'

  if (! all(domain %in% c("rr", "ri", "ir", "ii")) ) stop(err)

  if ( length(domain) == 1 ) {
    get(domain, x)
  } else if ( length(domain) == 2) {
    if ( use.cmplx1 ) {
      cmplx1(re = get(domain[1], x), im = get(domain[2], x))
    } else {
      complex(re = get(domain[1], x), im = get(domain[2], x))
    }
  } else if ( length(domain) == 3) {
    cmplx2(rr = get(domain[1], x), ri = get(domain[2], x), 
           ir = get(domain[3], x))
  } else if ( length(domain) == 4) {
    cmplx2(rr = get(domain[1], x), ri = get(domain[2], x), 
           ir = get(domain[3], x), ii = get(domain[4], x))
  } else {
    stop(err)
  }
}

#------------------------------------------------------------------------------
# Summary

#' @method as_tibble vctrs_cmplx2
#' @export
as_tibble.vctrs_cmplx2 <- function(x, ...) as_tibble(unclass(x))

#' @method summary vctrs_cmplx2
#' @export
summary.vctrs_cmplx2 <- function(object, ..., 
                                 digits = max(3, getOption("digits") - 3)) {
  summary(as_tibble(object))

}

#------------------------------------------------------------------------------
# Arithmetic

#---------------------------------------
# Boilerplate

#' @method vec_arith vctrs_cmplx2
#' @export
#' @export vec_arith.vctrs_cmplx2
vec_arith.vctrs_cmplx2 <- function(op, x, y, ...) {
  UseMethod("vec_arith.vctrs_cmplx2", y)
}

#' @method vec_arith.vctrs_cmplx2 default
#' @export
vec_arith.vctrs_cmplx2.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#---------------------------------------
# cmplx2

#' @method vec_arith.vctrs_cmplx2 vctrs_cmplx2
#' @export
vec_arith.vctrs_cmplx2.vctrs_cmplx2 <- function(op, x, y, ...) {
  switch(
    op,
    "+" = new_cmplx2(rr = x$rr + y$rr, ri = x$ri + y$ri,
                     ir = x$ir + y$ir, ii = x$ii + y$ii), 
    "-" = new_cmplx2(rr = x$rr - y$rr, ri = x$ri - y$ri,
                     ir = x$ir - y$ir, ii = x$ii - y$ii),
    stop_incompatible_op(op, x, y)
  )
}

#---------------------------------------
# complex 

#' @method vec_arith.vctrs_cmplx2 complex
#' @export
vec_arith.vctrs_cmplx2.complex <- function(op, x, y, ...) {
   vec_arith(op, x, vec_cast(y, new_cmplx2())) 
}

#' @method vec_arith.complex vctrs_cmplx2
#' @export
vec_arith.complex.vctrs_cmplx2 <- function(op, x, y, ...) {
  vec_arith(op, vec_cast(x, new_cmplx2()), y)
}

#---------------------------------------
# double

#' @method vec_arith.vctrs_cmplx2 double
#' @export
vec_arith.vctrs_cmplx2.double <- function(op, x, y, ...) {
  switch(
    op,
    "/" = new_cmplx2(rr = x$rr/y, ri = x$ri/y, ir = x$ir/y, ii = x$ii/y),
    vec_arith(op, x, cmplx2(rr = y))
  )
}

#' @method vec_arith.double vctrs_cmplx2
#' @export
vec_arith.double.vctrs_cmplx2 <- function(op, x, y, ...) {
  vec_arith(op, cmplx2(rr = x), y)
}



#------------------------------------------------------------------------------
# Math

#' @method vec_math vctrs_cmplx2
#' @export
vec_math.vctrs_cmplx2 <- function(.fn, .x, ...) {
  switch(.fn,
    sum = cmplx2(rr = sum(.x$rr), ri = sum(.x$ri), 
                 ir = sum(.x$ir), ii = sum(.x$ii)),
    mean = cmplx2(rr = sum(.x$rr), ri = sum(.x$ri), 
                  ir = sum(.x$ir), ii = sum(.x$ii))/as.numeric(length(.x)),
    vec_math_base(.fn, .x, ...)
  )
}

methods::setOldClass(c("vctrs_cmplx2", "vctrs_vctr"))
