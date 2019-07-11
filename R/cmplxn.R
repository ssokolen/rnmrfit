# Definition of a new set of complex classes that expand beyond the i notation


#' @import vctrs
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

#' @export
#' @method format vctrs_cmplx1
format.vctrs_cmplx1 <- function(x, ...) {
  r <- field(x, "r")
  i <- field(x, "i")
  
  signs <- c("1"="+", "0"="+", "-1"="")
  f <- function (x) signs[as.character(sign(x))]
  out <- paste0(r, f(i), i, "i")
  out[is.na(r) | is.na(i)] <- NA
  
  out
}

#' @export
#' @method vec_ptype_abbr vctrs_cmplx1
vec_ptype_abbr.vctrs_cmplx1 <- function(x) "cmplx1"

#' @export
#' @method vec_ptype_full vctrs_cmplx1
vec_ptype_full.vctrs_cmplx1 <- function(x) "complex1d"

#' @export
#' @method dim vctrs_cmplx1
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

#' @export
#' @method vec_ptype2 vctrs_cmplx1
vec_ptype2.vctrs_cmplx1 <- function(x, y, ...) {
  UseMethod("vec_ptype2.vctrs_cmplx1", y)
}

#' @export
#' @method vec_ptype2.vctrs_cmplx1 default
vec_ptype2.vctrs_cmplx1.default <- 
  function(x, y, ..., x_arg = "x", y_arg = "y") {
    vec_default_ptype2(x, y, x_arg = x_arg, y_arg = y_arg)
}

#' @export
#' @method vec_ptype2 complex
vec_ptype2.complex <- function(x, y, ...) {
  UseMethod("vec_ptype2.complex", y)
}

#' @export
#' @method vec_ptype2.complex default
vec_ptype2.complex.default <- 
  function(x, y, ..., x_arg = "x", y_arg = "y") {
    vec_default_ptype2(x, y, x_arg = x_arg, y_arg = y_arg)
}


#' @export
#' @method vec_ptype2.vctrs_cmplx1 vctrs_unspecified
vec_ptype2.vctrs_cmplx1.vctrs_unspecified <- function(x, y, ...) x

#' @export
#' @method vec_ptype2.vctrs_cmplx1 vctrs_cmplx1
vec_ptype2.vctrs_cmplx1.vctrs_cmplx1 <- function(x, y, ...) new_cmplx1()


#' @export
#' @method vec_ptype2.vctrs_cmplx1 complex
vec_ptype2.vctrs_cmplx1.complex <- function(x, y, ...) new_cmplx1()

#' @export
#' @method vec_ptype2.complex vctrs_cmplx1
vec_ptype2.complex.vctrs_cmplx1 <- function(x, y, ...) new_cmplx1()


#' @export
#' @method vec_ptype2.vctrs_cmplx1 double
vec_ptype2.vctrs_cmplx1.double <- function(x, y, ...) new_cmplx1()

#' @export
#' @method vec_ptype2.double vctrs_cmplx1
vec_ptype2.double.vctrs_cmplx1 <- function(x, y, ...) new_cmplx1()


#' @export
#' @method vec_ptype2.vctrs_cmplx1 data.frame
vec_ptype2.vctrs_cmplx1.data.frame<- function(x, y, ...) tibble()

#' @export
#' @method vec_ptype2.data.frame vctrs_cmplx1
vec_ptype2.data.frame.vctrs_cmplx1 <- function(x, y, ...) tibble()

#---------------------------------------
# Cast definitions

#' @export
#' @export vec_cast.vctrs_cmplx1
#' @method vec_cast vctrs_cmplx1
vec_cast.vctrs_cmplx1 <- function(x, to, ...) {
  UseMethod("vec_cast.vctrs_cmplx1")
}

#' @export
#' @method vec_cast.vctrs_cmplx1 default
vec_cast.vctrs_cmplx1.default <- function(x, to, ...) {
  vec_default_cast(x, to)
}


#' @export
#' @method vec_cast.vctrs_cmplx1 vctrs_cmplx1
vec_cast.vctrs_cmplx1.vctrs_cmplx1 <- function(x, to, ...) {
  x
}

#' @export
#' @method vec_cast.double vctrs_cmplx1
vec_cast.double.vctrs_cmplx1 <- function(x, to, ...) {
  Re(x)
}

#' @export
#' @method vec_cast.complex vctrs_cmplx1
vec_cast.complex.vctrs_cmplx1 <- function(x, to, ...) {
  complex(real = Re(x), imag = Im(x))
}

#' @export
#' @method vec_cast.data.frame vctrs_cmplx1
vec_cast.data.frame.vctrs_cmplx1 <- function(x, to, ...) {
  tibble(r = x$r, i = x$i)
}

#' @export
#' @method vec_cast.vctrs_cmplx1 double
vec_cast.vctrs_cmplx1.double <- function(x, to, ...) {
  cmplx1(r = x)
}

#' @export
#' @method vec_cast.vctrs_cmplx1 complex
vec_cast.vctrs_cmplx1.complex <- function(x, to, ...) {
  cmplx1(r = Re(x), i = Im(x))
}

#' @export
#' @method vec_cast.vctrs_cmplx1 data.frame
vec_cast.vctrs_cmplx1.data.frame <- function(x, to, ...) {
  cmplx1(r = x$r, i = x$i)
}

#---------------------------------------
# as_cmplx1 shortcut

#' @export
as_cmplx1 <- function(x) {
  vec_cast(x, new_cmplx1())
}

#---------------------------------------
#' Pack cmplx1 or cmplx2 column
#' 
#' This function undoes unpack() 
#'
#' @param x data.frame-like object.
#' @param column.name character specifying cmplx1 column in x.
#'
#' @return Modified x with new "r" and "i" columns.
#'
#' export
unpack <- function(x, column.name) {

  # x must have colnames
  err <- '"x" does not have colnames defined.'
  if ( identical(colnames(x), NULL) ) stop(err)

  # Check to make sure column exists and has vctrs_cmplx1 class
  err <- sprintf('"x" does not have a "%s" column', column.name)
  if (! column.name %in% colnames(x) ) stop(err) 

  err <- sprintf('"%s" column must be of type vctrs_cmplx1.', column.name)
  if (! 'vctrs_cmplx1' %in% class(x[[column.name]]) ) stop(err)

  index <- which(colnames(x) == column.name)
  column <- x[[column.name]]
  add_column(x[ , -index], r = column$r, i = column$i, .after = index)
}

#------------------------------------------------------------------------------
# Traditional Re()/Im()/Conj()

#' @export
#' @method Re vctrs_cmplx1
Re.vctrs_cmplx1 <- function(z) z$r

#' @export
#' @method Im vctrs_cmplx1
Im.vctrs_cmplx1 <- function(z) z$i

#' @export
#' @method Conj vctrs_cmplx1
Conj.vctrs_cmplx1 <- function(z) new_cmplx1(r = z$r, i = -z$i)

#------------------------------------------------------------------------------
# Summary

#' @export
#' @method as_tibble vctrs_cmplx1
as_tibble.vctrs_cmplx1 <- function(x, ...) as_tibble(unclass(x))

#' @export
#' @method summary vctrs_cmplx1
summary.vctrs_cmplx1 <- function(object, ..., 
                                 digits = max(3, getOption("digits") - 3)) {
  summary(as_tibble(object))

}

#------------------------------------------------------------------------------
# Arithmetic

#---------------------------------------
# Boilerplate

#' @export
#' @export vec_arith.vctrs_cmplx1
#' @method vec_arith vctrs_cmplx1
vec_arith.vctrs_cmplx1 <- function(op, x, y, ...) {
  UseMethod("vec_arith.vctrs_cmplx1", y)
}

#' @export
#' @method vec_arith.vctrs_cmplx1 default
vec_arith.vctrs_cmplx1.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#---------------------------------------
# cmplx1

#' @export
#' @method vec_arith.vctrs_cmplx1 vctrs_cmplx1
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
# double

#' @export
#' @method vec_arith.vctrs_cmplx1 numeric
vec_arith.vctrs_cmplx1.numeric <- function(op, x, y, ...) {
  switch(
    op,
    "/" = new_cmplx1(r = x$r/y, i = x$i/y),
    vec_arith(op, x, cmplx1(r = y))
  )
}

#' @export
#' @method vec_arith.numeric vctrs_cmplx1
vec_arith.numeric.vctrs_cmplx1 <- function(op, x, y, ...) {
  vec_arith(op, cmplx1(r = x), y)
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
vec_ptype2.vctrs_cmplx2.vctrs_unspecified <- function(x, y, ...) x

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
vec_ptype2.vctrs_cmplx2 <- function(x, y, ...) {
  UseMethod("vec_ptype2.vctrs_cmplx2", y)
}

vec_ptype2.vctrs_cmplx2.default <- function(x, y, ..., x_arg = "", y_arg = "") {
  stop_incompatible_type(x, y, x_arg = x_arg, y_arg = y_arg)
}



vec_ptype2.vctrs_cmplx2.vctrs_unspecified <- function(x, y, ...) x
vec_ptype2.vctrs_cmplx2.vctrs_cmplx2 <- function(x, y, ...) new_cmplx2()

vec_ptype2.vctrs_cmplx2.complex <- function(x, y, ...) new_cmplx2()
vec_ptype2.complex.vctrs_cmplx2 <- function(x, y, ...) new_cmplx2()

vec_ptype2.vctrs_cmplx2.double <- function(x, y, ...) new_cmplx2()
vec_ptype2.double.vctrs_cmplx2 <- function(x, y, ...) new_cmplx2()

#---------------------------------------
# Cast definitions

vec_cast.vctrs_cmplx2 <- function(x, to) UseMethod("vec_cast.vctrs_cmplx2")
vec_cast.vctrs_cmplx2.default <- function(x, to) vec_default_cast(x, to)

vec_cast.vctrs_cmplx2.vctrs_cmplx2 <- function(x, to) x
vec_cast.double.vctrs_cmplx2 <- function(x, to) Re(x)
vec_cast.complex.vctrs_cmplx2 <- function(x, to) {
  complex(real = Re(x), imag = Im(x))
}
vec_cast.vctrs_cmplx2.double <- function(x, to) cmplx2(rr = x)
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



#==============================================================================>
# Hack function to get around bind_rows() errors
#==============================================================================>



#---------------------------------------
#' Unpack/pack cmplx1 or cmplx2 column
#' 
#' Despite the conveniences offered by defining complex numbers as vctrs
#' objects, in some cases, it can still be more convenient to treat
#' real/imaginary values as separate values of a tibble data frame. One example
#' is when using the \code{group_by()} function from dplyr, as
#' \code{bind_rows()} has trouble with vctrs objects. \code{unpack()} and
#' \code{pack()} function provide a quick method to unpack/pack real/imginary
#' values into and out of separate columns.
#' 
#' @param x data.frame-like object containing cmplx1 or cmplx2 data.
#' @param column.name character specifying cmplx1 or cmplx2 column in x.
#' 
#' @return Modified tidyverse tibble. \code{unpack()} generates new "r" and "i"
#'         columns for cmplx1 data and "rr", "ri", "ir", "ii" columns for cmplx2
#'         data. \code{pack()} attempts to find the above columns and generates
#'         a new cmplx1 or cmplx2 column with the name specified by column.name.
#' 
#' @name unpack
#' @export
unpack <- function(x, column.name) {

  # x must have colnames
  err <- '"x" does not have colnames defined.'
  if ( identical(colnames(x), NULL) ) stop(err)

  # Check to make sure column exists and has vctrs_cmplx1 class
  err <- sprintf('"x" does not have a "%s" column', column.name)
  if (! column.name %in% colnames(x) ) stop(err) 

  err <- sprintf('"%s" column must be of type vctrs_cmplx1.', column.name)
  if (! 'vctrs_cmplx1' %in% class(x[[column.name]]) ) stop(err)

  index <- which(colnames(x) == column.name)

  wrn <- 'Multiple columns with "%s" name found, selecting the first.'
  if ( length(index) > 1 ) warning(wrn)
  index <- index[1]

  column <- x[[column.name]]
  new.columns <- as_tibble(column)

  if ( index == ncol(x) ) {
    unpacked <- cbind(x[, -ncol(x)], new.columns)
  } else if ( index == 1 ) {
    unpacked <- cbind(new.columns, x[, -1])
  } else {
    unpacked <- cbind(x[ , 1:(index-1)], new.columns, x[, index:ncol(x)])
  }

  new.index <- index:(index + ncol(new.columns) - 1)
  colnames(unpacked)[new.index] <- colnames(new.columns) 
  unpacked
}



#' @rdname unpack
#' @export
pack <- function(x, column.name) {

  x <- as_tibble(x)

  # Looking for appropriate columns
  cmplx1.names <- c('r', 'i')
  if ( all( cmplx1.names %in% colnames(x) ) ) {
    index <- which( colnames(x) %in% cmplx1.names )
    x[[index[1]]] <- as_cmplx1(x[, index])
    x <- x[, -index[-1]]
    colnames(x)[index[1]] <- column.name
  } else {
    err <- 'cmpplx1 or cmplx2 columns not found.'
    stop(err)
  }

  x
}

