# Definition of a super-class for 1D resonance data.



#==============================================================================>
#  NMRScaffold -- super-class for all 1D and 2D NMR classes
#==============================================================================>



#------------------------------------------------------------------------------
#' Super-class for all 1D and 2D peak descriptions.
#' 
#' This class is not meant to be used directly. Instead, it provides a set of
#' common methods for all 1D and 2D NMR objects. As such, it has no slots and
#' all of its methods are meant to be inherited by the objects listed above.
#' 
#' @name NMRScaffold-class
#' @export
NMRScaffold <- setClass("NMRScaffold", 
  contains = "VIRTUAL"
)



#==============================================================================>
#  Helper functions
#==============================================================================>



#---------------------------------------
#' Combine all components of object
#' 
#' This is an internal function used for all getter functions that output a
#' data.frame object -- stitching together the results of all components.
#' 
#' @param object An NMRScaffold object.
#' @param getter Getter function.
#' @param slot.name Slot name as character.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name combine_all
setGeneric("combine_all", 
  function(object, ...) standardGeneric("combine_all")
)

#' @rdname combine_all
setMethod("combine_all", "NMRScaffold", 
  function(object, getter, slot.name, ...) {

  slot.names <- slotNames(object)

  # The order here is important because NMRFit1D has both
  # lower.bounds/upper.bounds and children -- we want the children here

  # If the object has children pass the function on and combine
  if ( "children" %in% slot.names ) {
    out.list <- map(object@children, getter, include.id = TRUE)
    out <- do.call(rbind, out.list)
  }
  # If the object has corresponding slot, produce it
  else if ( slot.name %in% slot.names ) {
    out <- slot(object, slot.name)
  }
  # If the object has a dimensions slot, pass the function on and combine
  else if ( "dimensions" %in% slot.names ) {
    out <- combine_dimensions(object, getter, ...)
  } 
  # Otherwise, error
  else {
    err <- 'Object does not have "%s" slots or any valid components'
    err <- sprintf(err, slot.name)
    stop(err)
  }

  out
})

#---------------------------------------
#' Split input across all components of object
#' 
#' This is an internal function used for all setter functions that input a
#' data.frame object -- splitting apart input across all components.
#' 
#' @param object An NMRScaffold object.
#' @param setter Setter function.
#' @param slot.name Slot name as character.
#' @param value Value stored in specified slot.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name split_all
setGeneric("split_all", 
  function(object, setter, slot.name, value) standardGeneric("split_all")
)

#' @rdname split_all
#' @export
setMethod("split_all", "NMRScaffold", 
  function(object, setter, slot.name, value) {

  slot.names <- slotNames(object)

  # The setter is basically only used on peaks() so the order
  # issues with combine_all don't apply here

  # If the object has corresponding slot, assign as normal
  if ( slot.name %in% slot.names ) {
    slot(object, slot.name) <- value
    return(object)
  } 
  # If the object has a dimensions slot, pass the setter on
  else if ( "dimensions" %in% slot.names ) {
    return(split_dimensions(object, setter, value))
  } 

  # Otherwise, continue

  # First the input must be a data.frame of some sort
  err <- 'Input value must be a data.frame type object.'
  if (! 'data.frame' %in% class(value) ) stop(err)

  # Second the input must have child name column
  name <- object@children[[1]]@name
  err <- sprintf('Input data.frame must have a "%s" column.', name)
  if (! name %in% colnames(value) ) stop(err)

  # Third, the children column must only contain existing ids
  new.names <- unique(value[[name]])
  old.names <- unlist(lapply(object@children, id))
  logic <- new.names %in% old.names

  wrn <- sprintf('The following %s is not defined, ignoring: %s',
                 name, paste(new.names[!logic], collapse = ', '))
  if ( any(!logic) ) warning(wrn)

  # If all of the above is met, then split components
  f_select <- function(d) d[, ! grepl(name, colnames(d))]
  value.chunk <- by(value, value[[name]], f_select)
  indexes <- which(old.names %in% new.names)

  for ( i in indexes ) {
    child <- object@children[[i]]
    child <- setter(child, value.chunk[[old.names[i]]])
    object@children[[i]] <- child
  }
  
  object
})



#==============================================================================>
# Generics for basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
# Id

#---------------------------------------
#' Get object id
#' 
#' Generic convenience method to access the id of any NMRScaffold object.
#' 
#' @param object An NMRScaffold object.
#' @param ... Additional arguments passed to inheriting methods.
#'
#' @name id
#' @export
setGeneric("id", 
  function(object, ...) standardGeneric("id")
)

#' @rdname id
#' @export
setMethod("id", "NMRScaffold", 
  function(object) object@id
)



#------------------------------------------------------------------------------
# Peaks

#---------------------------------------
#' Get object peaks
#' 
#' Generic convenience method to access the peak definitions of any
#' NMRScaffod object. If used on anything more complicated than an
#' NMRResonance1D, the function combines the peaks data frames of component
#' resonances.
#' 
#' @param object An NMRScaffold object.
#' @param include.id TRUE to return a column of component ids.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name peaks
#' @export
setGeneric("peaks", 
  function(object, ...) standardGeneric("peaks")
)

#' @rdname peaks
#' @export
setMethod("peaks", "NMRScaffold", 
  function(object, include.id = FALSE) {
  
  out <- combine_all(object, peaks, "peaks")
  if ( include.id & (nrow(out) > 0) ) {
    out <- cbind(id = object@id, out)
    colnames(out)[1] <- object@name
  }
  out

})

#---------------------------------------
#' Set object peaks
#' 
#' Generic convenience method to set the peak definitions of any NMRScaffold
#' object. This is primarily intended as an internal method, so use with
#' caution. Changing peak definitions after an object has been defined may have
#' unpredictable consequences.
#' 
#' @param object An NMRScaffold object.
#' @param value A data frame with a minimum of "position", "width", "height",
#'              and "fraction.gauss" columns. Peaks may be defined by one to
#'              four columns of "peak", "resonance", "species", and "dimension"
#'              depending on the nature of the original object.
#' 
#' @name peaks-set
#' @export
setGeneric("peaks<-", 
  function(object, value) standardGeneric("peaks<-")
)

#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRScaffold",
  function(object, value) split_all(object, `peaks<-`, "peaks", value)
)


#------------------------------------------------------------------------------
# Couplings

#---------------------------------------
#' Get object couplings 
#' 
#' Generic convenience method to access the coupling definitions of any
#' NMRScaffold object. If used on anything more complicated than an
#' NMRResonance1D, the function combines the couplings data frames of component
#' resonances.
#' 
#' @param object An NMRScaffold object.
#' @param include.id TRUE to return a column of resonance or species ids.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name couplings
#' @export
setGeneric("couplings", 
  function(object, ...) standardGeneric("couplings")
)

#' @rdname couplings
#' @export
setMethod("couplings", "NMRScaffold", 
  function(object, include.id = FALSE) {

  out <- combine_all(object, couplings, "couplings")
  if ( include.id & (nrow(out) > 0) ) {
    out <- cbind(id1 = object@id, id2 = object@id, out)
    colnames(out)[1:2] <- paste(object@name, 1:2, sep = '.')
  }
  out

})



#------------------------------------------------------------------------------
# Bounds

#---------------------------------------
#' Get object bounds
#' 
#' Generic convenience method to access the bounds of any NMRScaffold object.
#' If used on anything more complicated than an NMRResonance1D, the function
#' combines the bounds data frames of component resonances. Whereas
#' lower_bounds() and upper_bounds() return a data.frame to match peaks(),
#' bounds() groups the lower and upper bounds into a list.
#' 
#' @param object An NMRScaffold object.
#' @param include.id TRUE to return a column of resonance or species ids.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name bounds
#' @export
setGeneric("bounds", 
  function(object, ...) standardGeneric("bounds")
)

#' @rdname bounds
#' @export
setGeneric("lower_bounds", 
  function(object, ...) standardGeneric("lower_bounds")
)

#' @rdname bounds
#' @export
setGeneric("upper_bounds", 
  function(object, ...) standardGeneric("upper_bounds")
)

#' @rdname bounds 
#' @export
setMethod("bounds", "NMRScaffold", 
  function(object, include.id = FALSE) {
  
    list(lower = lower_bounds(object, include.id),
         upper = upper_bounds(object, include.id))

})

#' @rdname bounds
#' @export
setMethod("lower_bounds", "NMRScaffold", 
  function(object, include.id = FALSE) {

  out <- combine_all(object, lower_bounds, "lower.bounds")
  if ( include.id & (nrow(out) > 0) ) {
    out <- cbind(id = object@id, out)
    colnames(out)[1] <- object@name
  }
  out

})

#' @rdname bounds
#' @export
setMethod("upper_bounds", "NMRScaffold", 
  function(object, include.id = FALSE) {

  out <- combine_all(object, upper_bounds, "upper.bounds")
  if ( include.id & (nrow(out) > 0) ) {
    out <- cbind(id = object@id, out)
    colnames(out)[1] <- object@name
  }
  out

})



#==============================================================================>
# Validation methods 
#==============================================================================>



#---------------------------------------
#' Check data conformity
#' 
#' This is intended to a primarily internal function used to check whether
#' provided data conforms to the scaffold. The check is currently limited to
#' ensuring the correct number of dimensions and that data covers the chemical
#' shift of specific positions
#' 
#' @param object An NMRScaffold object.
#' @param nmrdata An NMRData object.
#' @param error True to issue an error message, otherwise, return message as
#'              string.
#' @param ... Additional arguments passed to inheriting methods.
#'
#' @return TRUE if data conforms, FALSE otherwise.
#' 
#' @name check_conformity
#' @export
setGeneric("check_conformity", 
  function(object, ...) standardGeneric("check_conformity")
)



#==============================================================================>
# Convenience methods for bounds
#==============================================================================>



#---------------------------------------
#' Update bounds
#' 
#' This is an internal function used to update all the bounds of NMRScaffold object at once while ensuring that they are not widened if that is not desired. Do not use this function directly.
#' 
#' @param object An NMRScaffold object.
#' @param lower.bounds data.frame of lower bounds.
#' @param upper.bounds data.frame of upper bounds.
#' @param widen FALSE to prevent new bounds from widening existing bounds.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name update_bounds
setGeneric("update_bounds", 
  function(object, ...) standardGeneric("update_bounds")
)

#' @rdname update_bounds
setMethod("update_bounds", "NMRScaffold", 
  function(object, lower.bounds, upper.bounds, widen) {
    
  # Ensuring that parameters are only widened if desired
  columns <- c('position', 'width', 'height', 'fraction.gauss')

  new.lower <- unlist(lower.bounds[ , columns])
  old.lower <- unlist(lower_bounds(object)[ , columns])

  new.upper <- unlist(upper.bounds[ , columns])
  old.upper <- unlist(upper_bounds(object)[ , columns])
  
  if (! widen ) {
    new.lower <- ifelse(new.lower < old.lower, old.lower, new.lower)
    new.upper <- ifelse(new.upper > old.upper, old.upper, new.upper)
  }

  object@lower.bounds[ , columns] <- new.lower
  object@upper.bounds[ , columns] <- new.upper

  object
})



#------------------------------------------------------------------------------
#' Set general bounds of an NMRScaffold object
#' 
#' This function provides a convenience method for generating bounds using a
#' simple set of lower and upper constraints on basic peak parameters such as
#' position, height, width, and fraction.gauss. The term "general" refers to
#' the fact that the same bounds are applied to each and every peak, regardless
#' of current parameter values. These bounds can be normalized to a set of data
#' using the optional nmrdata argument. If used on an anything more complicated
#' than an NMRResonance1D, the function propogates the bounds to component
#' resonances.
#' 
#' In practice, general bounds are primarily useful for placing a hard
#' constraint on peak widths and preventing negative heights. Values of 0 for
#' widths and sometimes height can also cause issues during optimization, so
#' simple general bounds can be used to prevent errors.
#' 
#' @param object An NMRScaffold object.
#' @param position A vector of two elements corresponding to a lower and upper
#'                 bound for peak position. If nmrdata is provided, 0
#'                 corresponds to the leftmost range of the data and 1 to the
#'                 rightmost. Otherwise, the units are in ppm.
#' @param height A vector of two elements corresponding to a lower and upper
#'               bound for peak height. If nmrdata is provided, 0 corresponds to
#'               the lowest value of spectral intensity (in the real domain) and
#'               1 to the largest. Otherwise, the units correspond to arbitrary
#'               spectral intensity values.
#' @param width A vector of two elements corresponding to a lower and upper
#'              bound for peak width in Hz. If nmrdata is provided, values are
#'              taken as fraction of the general data range. So 0.1 would
#'              correspond to a nominal peak width that covers a tenth of the
#'              general data range.
#' @param fraction.gauss A vector of two elements corresponding to a lower and
#'                       upper bound for the Gaussian fraction of the peak. This
#'                       can be set to c(0, 0) to force Lorentzian peaks. Any
#'                       values smaller than 0 will be treated as 0 and any
#'                       values greater than 1 will be treated as 1.
#' @param baseline Only applicable to NMRFit1D objects. A complex vector of two
#'                 elements corresponding to a lower and upper bound for both
#'                 real and imaginary baseline control points (which roughly
#'                 corresponds to a baseline value). If the vector has no
#'                 imaginary component, the imaginary baseline is left
#'                 unbounded. If nmrdata is provided, 0 corresponds to the
#'                 lowest value of spectral intensity and 1 to the largest.
#'                 Otherwise, the units correspond to arbitrary spectral
#'                 intensity values.
#' @param phase Only applicable to NMRFit1D objects. A vector of two elements
#'              corresponding to a lower and upper bound for the phase
#'              correction (in radians) at any point in the chemical shift
#'              range.
#' @param nmrdata An optional NMRData1D object that can serve as a reference
#'                point for the bounds.
#' @param widen FALSE to prevent new bounds from widening existing bounds.
#' @param dimensions Bounds are projected to all dimensions by default. Specific
#'                   dimensions can be specified as "direct" or "indirect".
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A new object with modified bounds.
#' 
#' @name set_general_bounds
#' @export
setGeneric("set_general_bounds", 
  function(object, ...) standardGeneric("set_general_bounds")
)

#' @rdname set_general_bounds
#' @export
setMethod("set_general_bounds", "NMRScaffold",
  function(object, position = NULL, height = NULL, width = NULL, 
           fraction.gauss = NULL, nmrdata = NULL, widen = FALSE,
           dimensions = NULL) {
    
  # If object has children, just apply function recursively
  if ( "children" %in% slotNames(object) ) {

    object@children <- lapply(object@children, set_general_bounds, 
      position = position, height = height, width = width, 
      fraction.gauss = fraction.gauss, nmrdata = nmrdata, widen = widen
    )
    return(object)
  }

  # If object has dimensions, propagate as necessary
  if ( "dimensions" %in% slotNames(object) ) {

    valid.dimensions <- c("direct", "indirect")
    if ( is.null(dimensions) ) dimensions <- valid.dimensions

    if (! all(dimensions %in% valid.dimensions) ) {
      err <- paste('"dimensions" must be one of', 
                   paste(dimensions, collapse = ", "))
      stop(err)
    }

    f_bounds <- function(child, dimension) {
      if (! is.null(nmrdata) ) {
        nmrdata <- projection(nmrdata, dimension)
      }
      set_general_bounds(child, position, height, width,
                         fraction.gauss, nmrdata, widen)
    }
    
    for ( name in dimensions ) {
      object@dimensions[[name]] <- f_bounds(object@dimensions[[name]], name)
    }

    return(object)
  }

  # Otherwise, continue

  lower <- object@lower.bounds
  upper <- object@upper.bounds

  #---------------------------------------
  # Scaling all bounds if nmrdata has been provided
  if (! is.null(nmrdata) ) {

    if ( class(nmrdata) != 'NMRData1D' ) {
      err <- '"nmrdata" must be a valid NMRData1D object.'
      stop(err)
    }

    processed <- nmrdata@processed
    y.range <- max(Re(processed$intensity)) - min(Re(processed$intensity))
    x.range <- max(processed$direct.shift) - min(processed$direct.shift)

    position <- position * x.range + min(processed$direct.shift)
    height <- height * y.range

    sfo1 <- get_parameter(nmrdata, 'sfo1', 'acqus', error = TRUE)
    width <- width * (x.range * sfo1)
  }

  #---------------------------------------
  # Defining a bound check function
  .check_bounds <- function(bounds) {

    if ( length(bounds) != 2 ) {
      err <- paste("All bounds must be vectors of two elements consisting",
                   "of a lower and upper bound.")
      stop(err)
    }

    if ( bounds[1] > bounds[2] ) {
      err <- "Lower bound must be smaller than upper bound."
      stop(err)
    }

  }

  #---------------------------------------
  # Creating a list of bounds to loop through each in term
  bounds = list(position = position, height = height, width = width,
                fraction.gauss = fraction.gauss)

  for ( parameter in names(bounds) ) {
    if ( ! is.null(bounds[[parameter]]) ) {
      .check_bounds(bounds[[parameter]])
      lower[[parameter]] <- bounds[[parameter]][1]
      upper[[parameter]] <- bounds[[parameter]][2]
    }
  }

  # Fraction gauss is a little different because it must be 0-1
  lower$fraction.gauss[lower$fraction.gauss < 0] <- 0
  upper$fraction.gauss[upper$fraction.gauss > 1] <- 1

  # Ensuring that parameters are only widened if desired
  update_bounds(object, lower, upper, widen)
})


#------------------------------------------------------------------------------
#' Set offset bounds of an NMRScaffold object
#' 
#' This function provides a convenience method for generating bounds using a
#' simple set of lower and upper constraints on basic peak parameters such as
#' position, height, width, fraction.gauss. The term "offset" refers to the
#' fact that bounds are applied as an offset to current values of the
#' parameters. These bounds can be expressed in absolute (e.g. -0.1 and +0.1
#' ppm) or relative (e.g. -1 percent and +1 percent) terms. If used on an
#' anything more complicated than an NMRResonance1D, the function propogates
#' the bounds to component resonances.
#' 
#' In practice, offset bounds are primarily useful for preventing peak
#' positions from drifting too much from initial guesses and for fine-tuning a
#' fit once an initial optimization is performed. It is not recommended to use
#' strict offset bounds based on rough initial parameter guesses.
#' 
#' @param object An NMRScaffold object.
#' @param position A vector of two elements to be added to current peak
#'                 positions to generate a set of lower and upper bounds. If
#'                 relative is true, the values are treated as fractions to be
#'                 multipled by the current peak position before addition.
#'                 fraction of the current position.
#' @param height A vector of two elements to be added to current peak heights to
#'               generate a set of lower and upper bounds. If relative is true,
#'               the values are treated as fractions to be multipled by the
#'               current peak heights before addition.
#' @param width A vector of two elements to be added to current peak widths to
#'              generate a set of lower and upper bounds. If relative is true,
#'              the values are treated as fractions to be multipled by the
#'              current peak heights before addition.
#' @param relative TRUE to treat values as relative fractions, FALSE to apply
#'                 them directly.
#' @param widen FALSE to prevent new bounds from widening existing bounds.
#' @param dimensions Bounds are projected to all dimensions by default. Specific
#'                   dimensions can be specified as "direct" or "indirect".
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A new object with modified bounds.
#' 
#' @name set_offset_bounds
#' @export
setGeneric("set_offset_bounds", 
  function(object, ...) standardGeneric("set_offset_bounds")
)

#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRScaffold",
  function(object, position = NULL, height = NULL, width = NULL, 
           relative = FALSE, widen = FALSE, dimensions = NULL) {
    
  # If object has children, just apply function recursively
  if ( "children" %in% slotNames(object) ) {
    object@children <- lapply(object@children, set_offset_bounds, 
      position, height, width, relative, widen
    )
    return(object)
  }

  # If object has dimensions, propagate as necessary
  if ( "dimensions" %in% slotNames(object) ) {

    valid.dimensions <- c("direct", "indirect")
    if ( is.null(dimensions) ) dimensions <- valid.dimensions

    if (! all(dimensions %in% valid.dimensions) ) {
      err <- paste('"dimensions" must be one of', 
                   paste(dimensions, collapse = ", "))
      stop(err)
    }

    f_bounds <- function(child, dimension) {
      set_offset_bounds(child, position, height, width,
                        relative, widen)
    }
    
    for ( name in dimensions ) {
      object@dimensions[[name]] <- f_bounds(object@dimensions[[name]], name)
    }

    return(object)
  }


  # Otherwise, continue

  peaks <- object@peaks
  lower <- object@lower.bounds
  upper <- object@upper.bounds

  #---------------------------------------
  # Defining a bound check function
  .check_bounds <- function(bounds) {

    if ( length(bounds) != 2 ) {
      err <- paste("All bounds must be vectors of two elements consisting",
                   "of a lower and upper bound.")
      stop(err)
    }

    if ( bounds[1] > 0 ) {
      err <- paste("Lower offsets must be negative so that resulting bounds",
                   "include initial values.")
      stop(err)
    }

    if ( bounds[2] < 0 ) {
      err <- paste("Upper offsets must be positive so that resulting bounds",
                   "include initial values.")
      stop(err)
    }

    if ( bounds[1] > bounds[2] ) {
      err <- "Lower bound must be smaller than upper bound."
      stop(err)
    }

  }

  #---------------------------------------
  # Creating a list of bounds to loop through each in term
  bounds = list(position = position, height = height, width = width)

  for ( parameter in names(bounds) ) {
    if ( length(bounds[[parameter]]) > 0 ) {
      .check_bounds(bounds[[parameter]])

      lower.offset <- bounds[[parameter]][1]
      upper.offset <- bounds[[parameter]][2]
      
     if ( relative ) {
        lower.offset <- lower.offset*peaks[[parameter]]
        upper.offset <- upper.offset*peaks[[parameter]]
      } 

      lower[[parameter]] <- peaks[[parameter]] + lower.offset
      upper[[parameter]] <- peaks[[parameter]] + upper.offset
    }
  }

  # Ensuring that parameters are only widened if desired
  update_bounds(object, lower, upper, widen)
})



#------------------------------------------------------------------------------
#' Set conservative bounds on an NMRScaffold object
#' 
#' A convenience function that sets reasonable bounds on an NMRScaffold
#' object. These bounds are assumed to be widely applicable to most simple NMR
#' data. Each set of bounds can be turned on or off as necessary. A slightly
#' better set of bounds can be selected if a reference NMRData1D object is
#' provided. If used on an anything more complicated than an NMRResonance1D,
#' the function propogates the bounds to component resonances.
#' 
#' @param object An NMRScaffold object.
#' @param position Without reference data, position is limited to plus or minus
#'                 0.1 ppm. With reference data, the position of the peaks is
#'                 also forced inside the domain of the data. FALSE to disable.
#' @param height Without reference data, height is set to strictly positive.
#'               With reference data, height is also limited to no more than 150
#'               percent of the maximum peak value. FALSE to disable.
#' @param width Without reference data, minimum peak width is set to almost, but
#'              not quite 0 Hz (1e-3 Hz) and a maximum peak width of 3 Hz. With
#'              reference data, peak width is also prevented from being more
#'              than 20 percent of the data range. FALSE to disable.
#' @param baseline Without reference data, the baseline control points are left
#'                 unrestricted. With reference data, the real baseline control
#'                 points are limited to no more than 50 percent of the maximum
#'                 peak value. FALSE to disable.
#' @param phase With or without reference data, phase correction is limited to
#'              -pi/2 to pi/2.
#' @param nmrdata An optional NMRData1D object that can serve as a reference
#'                point for the bounds.
#' @param widen FALSE to prevent new bounds from widening existing bounds.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A new object with modified bounds.
#' 
#' @name set_conservative_bounds
#' @export
setGeneric("set_conservative_bounds", 
  function(object, ...) standardGeneric("set_conservative_bounds")
)

#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRScaffold",
  function(object, position = TRUE,  height = TRUE, width = TRUE, 
           nmrdata = NULL, widen = FALSE) { 

  # First, do a single pass over general bounds with no reference
  if ( height )  gen.height <- c(0, Inf)
  else gen.height <- NULL

  if ( width ) gen.width <- c(0.003, 5)
  else gen.width <- NULL

  object <- set_general_bounds(object, height = gen.height, width = gen.width,
                               widen = widen)

  # Adding position offsets
  if ( position ) {
    object <- set_offset_bounds(object, position = c(-0.1, 0.1), widen = widen)
  }

  # If nmrdata is provided, add further constraints  
  if (! is.null(nmrdata) ) {
    
    if (! class(nmrdata) %in% c('NMRData1D', 'NMRData2D') ) {
      err <- '"nmrdata" must be a valid NMRData1D object.'
      stop(err)
    } 

    if ( position )  gen.position <- c(0, 1)
    else gen.position <- NULL

    if ( height ) gen.height <- c(0, 1.5)
    else gen.height <- NULL

    if ( width ) gen.width <- c(0, 0.2)
    else gen.width <- NULL

    object <- set_general_bounds(object, position = gen.position, 
                                 height = gen.height, width = gen.width,
                                 nmrdata = nmrdata, widen = widen)
  }

  object
})



#------------------------------------------------------------------------------
#' Set peak type of an NMRScaffold object
#' 
#' The peak type of an NMRScaffold object is governed by the fraction.gauss
#' parameter and can fluidly go from pure Lorentz to Gauss via the combined
#' Voigt lineshape. However, it can be computationally efficient to force a
#' specific peak type. This function provides a shortcut for doing so by
#' changing the fraction.gauss parameter as well as lower and upper bounds.
#' 
#' @param object An NMRScaffold object.
#' @param peak.type One of either "lorentz", "voigt", "gauss", or "any" where
#'                  "any" clears existing bounds on the fraction.gauss
#'                  parameter.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A new object with modified bounds and peak parameters.
#' 
#' @name set_peak_type
#' @export
setGeneric("set_peak_type", 
  function(object, ...) standardGeneric("set_peak_type")
)

#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRScaffold",
  function(object, peak.type) {
    
  # If object has children, just apply function recursively
  if ( "children" %in% slotNames(object) ) {
    object@children <- lapply(object@children, set_peak_type, peak.type) 
    return(object)
  }

  # If object has dimensions, propagate as necessary
  if ( "dimensions" %in% slotNames(object) ) {

    valid.dimensions <- c("direct", "indirect")
    if ( is.null(dimensions) ) dimensions <- valid.dimensions

    if (! all(dimensions %in% valid.dimensions) ) {
      err <- paste('"dimensions" must be one of', 
                   paste(dimensions, collapse = ", "))
      stop(err)
    }

    f_bounds <- function(child, dimension) {
      set_peak_type(child, peak.type)
    }
    
    for ( name in dimensions ) {
      object@dimensions[[name]] <- f_bounds(object@dimensions[[name]], name)
    }

    return(object)
  }

  # Otherwise, continue

  peaks <- object@peaks
  lower <- object@lower.bounds
  upper <- object@upper.bounds

  # Getting rid of empty spaces and capitals
  peak.type <- tolower(gsub('\\s', '', peak.type))
  peak.types <- c('lorentz', 'voigt', 'gauss', 'any')
  peak.type <- pmatch(peak.type, peak.types, nomatch = -1)

  if ( peak.type == 1 ) {
    lower$fraction.gauss <- 0
    upper$fraction.gauss <- 0
    peaks$fraction.gauss <- 0
  } else if ( peak.type == 2 ) {
    lower$fraction.gauss <- 1e-6
    upper$fraction.gauss <- 1 - 1e-6
    peaks$fraction.gauss <- 0.5
  } else if ( peak.type == 3 ) {
    lower$fraction.gauss <- 1
    upper$fraction.gauss <- 1
    peaks$fraction.gauss <- 1
  } else if ( peak.type == 4 ) {
    lower$fraction.gauss <- 0
    upper$fraction.gauss <- 1
  } else {
    peak.types <- paste(peak.types, collapse = ', ')
    err <- sprintf('Peak type must be one of %s', peak.types)
    stop(err)
  }

  object@lower.bounds <- lower
  object@upper.bounds <- upper
  object@peaks <- peaks

  object
})



#==============================================================================>
# Fit Methods
#==============================================================================>



#------------------------------------------------------------------------------
#' Parse constraints for fitting
#' 
#' This is an internal function used to parse equality and inequality
#' constraints into a form that can be used for lower level fit code.
#' Essentially, each constraint is defined as a vector with a variable number
#' of elements codified as follows: 1) Code of either 0 (position), 1 (width),
#' 2 (height), 3 (fraction.gauss), 4 (area); 2) The equality or inequality
#' target value; 3+) Integers corresponding to peak indexes in the overall
#' parameter vector. Positive vs negative values are treated differently
#' depending on the code. For positions, parameters with negative indexes are
#' subtracted while for width and area, all parameters with negative indexes
#' are added up before dividing the positive ones.
#' 
#' @param object An NMRScaffold object.
#' @param direct.span The ppm range of the fit used to scale position
#'                    differences in the direct dimension.
#' @param indirect.span The ppm range of the fit used to scale position
#'                      differences in the indirect dimension.
#' 
#' @name parse_constraints
setGeneric("parse_constraints", 
  function(object, ...) {
    standardGeneric("parse_constraints")
})

#' @rdname parse_constraints
#' @export
setMethod("parse_constraints", "NMRScaffold",
  function(object, direct.span = 1, indirect.span = 1) {

    #---------------------------------------
    # Defining separate functions for dealing with the 5 different constraints.
    
    # Differences in peak position
    f_position <- function(leeway, difference, indexes) {
      if ( leeway == 0 ) {
        list(c(0, difference, indexes))
      } else {
        list(c(0, difference*(1+leeway), indexes),
             c(0, -difference*(1-leeway), -indexes))
      }
    }

    # height currently ignored

    # Ratios of peak width (fixed at equal)
    f_width <- function(leeway, indexes) {
      if ( leeway == 0 ) {
        list(c(1, 1, indexes))
      } else {
        list(c(1, (1+leeway), indexes),
             c(1, 1/(1-leeway), -indexes))
      }
    }

    # Differences of fraction.gauss (fixed at equal)
    f_fraction <- function(leeway, indexes) {
      if ( leeway == 0 ) {
        list(c(3, 0, indexes))
      } else {
        list(c(3, leeway, indexes),
             c(3, leeway, -indexes))
      }
    }

    # Ratios of area
    f_area <- function(leeway, ratio, indexes) {
      if ( leeway == 0 ) {
        list(c(4, ratio, indexes))
      } else {
        # All ratios are converted to be greater than 1 to allow for
        # standardized lower and upper bounds
        if ( ratio < 1 ) {
          ratio <- 1/ratio
          indexes <- -indexes
        }

        list(c(4, ratio*(1+leeway),  indexes),
             c(4, 1/(ratio*(1-leeway)), -indexes))
      }
    }

    #---------------------------------------
    # Applying the constraints

    eq.constraints <- list()
    ineq.constraints <- list()

    #---------------------------------------
    # Normalizing data
    # (technically, the following should only be applied to mixture
    # objects, but it's pretty straightforward to generate a trivial mixture
    # from resonance/species -- not tested yet, but should work)

    if ( class(object) %in% c("NMRResonance1D", "NMRSpecies1D") ) {
      object <- nmrmixture_1d(object)
    } else if ( class(object) %in% c("NMRResonance2D", "NMRSpecies2D") ) {
      object <- nmrmixture_2d(object)
    }
    
    peaks <- peaks(object)
    peaks$resonance <- as.character(peaks$resonance)

    # If we are dealing with 1D data, add placeholder dimension column
    if (! "dimension" %in% colnames(peaks) ) {
      peaks$dimensions = "direct"
    }
    
    # Combine direct and indirect span
    x.span <- list(direct = direct.span, indirect = indirect.span)

    #---------------------------------------
    # Generating constraints

    # To ensure that individual leeway values are considered, looping through
    # each species/resonance one at a time
    for ( specie in object@children ) {

      # At the species level, constraints are based on overall area sums
      # (only handled in the direct dimensions with the indirect dimension
      # being automatically handled by a more general constraint on height)
      leeway <- abs(specie@connections.leeway)

      # First dealing with conenctions between resonances
      connections <- specie@connections
      
      if ( nrow(connections) > 0 ) {
        
        for ( i in 1:nrow(connections) ) {

          resonance.1 <- connections$resonance.1[i]
          resonance.2 <- connections$resonance.2[i]
          ratio <- connections$area.ratio[i]

          logic <- specie@id == peaks$species &
                   (peaks$dimension == "direct")
          logic.1 <- logic & (resonance.1 == peaks$resonance)
          logic.2 <- logic & (resonance.2 == peaks$resonance)
          indexes <- c(which(logic.2), -which(logic.1))

          if ( leeway == 0 ) {
            eq.constraints <- c(eq.constraints, 
                                f_area(leeway, ratio, indexes))
          } else {
            ineq.constraints <- c(ineq.constraints, 
                                f_area(leeway, ratio, indexes))
          }
          
        }
      }

      # Then looping through each specific resonance within each species
      for ( resonance in specie@children ) {

        # If dealing with 2D data, use dimensions list, otherwise fake one
        if ( "dimensions" %in% slotNames(resonance) ) {
          dimensions <- resonance@dimensions
        } else {
          dimensions <- list(direct = resonance)
        }

        # Loop over dimensions
        for ( dimension in names(dimensions) ) {

          # At the species level, constraints are based on overall area sums
          id <- dimensions[[dimension]]@id
          couplings.leeway <- dimensions[[dimension]]@couplings.leeway

          position.leeway <- abs(couplings.leeway$position)
          width.leeway <- abs(couplings.leeway$width)
          fraction.leeway <- abs(couplings.leeway$fraction.gauss)
          area.leeway <- abs(couplings.leeway$area)

          # Only continue if there are couplings defined
          couplings <- dimensions[[dimension]]@couplings

          if ( nrow(couplings) > 0 ) {
            
            for ( i in 1:nrow(couplings) ) {
              peak.1 <- couplings$peak.1[i]
              peak.2 <- couplings$peak.2[i]
              diff <- couplings$position.difference[i]/x.span[[dimension]]
              ratio <- couplings$area.ratio[i]

              logic <- (specie@id == peaks$species) & 
                       (resonance@id == peaks$resonance) &
                       (dimension == peaks$dimension)

              logic.1 <- logic & (peak.1 == peaks$peak)
              logic.2 <- logic & (peak.2 == peaks$peak)
              indexes <- c(which(logic.2), -which(logic.1))

              # Each coupling constraint includes position, width, and area

              # First, the position
              leeway <- position.leeway
              if ( leeway == 0 ) {
                eq.constraints <- c(eq.constraints, 
                                    f_position(leeway, diff, indexes))
              } else {
                ineq.constraints <- c(ineq.constraints, 
                                    f_position(leeway, diff, indexes))
              }

              # Then, the width (same by default)
              leeway <- width.leeway
              if ( leeway == 0 ) {
                eq.constraints <- c(eq.constraints, 
                                    f_width(leeway, indexes))
              } else {
                ineq.constraints <- c(ineq.constraints, 
                                    f_width(leeway, indexes))
              }

              # Then, the fraction Gauss (same by default)
              leeway <- fraction.leeway
              if ( leeway == 0 ) {
                eq.constraints <- c(eq.constraints, 
                                    f_fraction(leeway, indexes))
              } else {
                ineq.constraints <- c(ineq.constraints, 
                                    f_fraction(leeway, indexes))
              }

              # Finally, the area
              leeway <- area.leeway
              if ( leeway == 0 ) {
                eq.constraints <- c(eq.constraints, 
                                    f_area(leeway, ratio, indexes))
              } else {
                ineq.constraints <- c(ineq.constraints, 
                                    f_area(leeway, ratio, indexes))
              }
            }
          }
        }
      }
    }

    # Returning list of constraints
    list(eq.constraints, ineq.constraints)
})
