# Definition of a super-class for 1D resonance data.



#==============================================================================>
#  NMRScaffold1D -- super-class for Resonance, Species, Mixture and Fit
#==============================================================================>



#------------------------------------------------------------------------------
#' Super-class for all 1D peak descriptions.
#' 
#' This class is not meant to be used directly. Instead, it provides a common
#' framework for methods around visualization and inspection of NMRResonance1D,
#' NMRSpecies1D, NMRMixture1D, and NMRFit1D objects. As such, it has no slots
#' and all of its methods are meant to be inherited by the objects listed
#' above.
#' 
#' @name NMRScaffold1D-class
#' @export
NMRScaffold1D <- setClass("NMRScaffold1D", 
  contains = "VIRTUAL"
)



#==============================================================================>
#  Display function
#==============================================================================>



#------------------------------------------------------------------------------
#' Display NMRScaffold1D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRScaffold1D", 
  function(object) {

  # Generating compiled data frames
  id <- ifelse( 'id' %in% slotNames(object), sprintf('(%s)', id(object)), '')
  peaks <- peaks(object)
  lower <- lower_bounds(object)
  upper <- upper_bounds(object)
  couplings <- couplings(object)

  # Id (if it exists)
  cat(sprintf('An object of %s class %s\n\n', class(object), id))

  # Peaks
  cat('Peaks:\n\n')
  print(peaks)
  cat('\n')

  # Bounds
  columns <- c('position', 'width', 'height', 'fraction.gauss')

  range <- paste('(', lower, ', ', upper, ')', sep = '')
  peaks[ , columns] <- range

  cat('Bounds (lower, upper):\n\n')
  print(peaks)
  cat('\n')   

  # Couplings
  if ( nrow(couplings) > 0 ) {
    cat('Couplings:\n\n')
    print(couplings)
    cat('\n')
  }
  else {
    cat('No couplings defined.\n')
  }

})



#==============================================================================>
#  Helper functions
#==============================================================================>



#---------------------------------------
#' Combine children if they exist
#' 
#' This is an internal function used for all getter functions that output a
#' data.frame object. If the object is a species or mixture that has children,
#' the getter is applied to all children, a children column is added and the
#' resulting objects are stitched together.
#' 
#' @param object An NMRScaffold1D object.
#' @param getter Getter function.
#' @param slot.name Slot name as character.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name combine_childrens
setGeneric("combine_children", 
  function(object, ...) standardGeneric("combine_children")
)

#' @rdname combine_children
setMethod("combine_children", "NMRScaffold1D", 
  function(object, getter, slot.name, ...) {
  
  if ( "children" %in% slotNames(object) ) {
    out.list <- lapply(object@children, getter, include.id = TRUE)
    out <- do.call(rbind, out.list)
  }
  else if ( slot.name %in% slotNames(object) ) {
    out <- slot(object, slot.name)
  } else {
    err <- 'Object does not have "%s" or "children" slots'
    err <- sprintf(err, slot.name)
    stop(err)
  }

  out
})

#---------------------------------------
#' Split children if they exist
#' 
#' This is an internal function used for all setter functions that input a
#' data.frame object. Essentially, the input value is first split based on a
#' "resonance" or "species" column, with the setter being passed on to the 
#' specified children.
#' 
#' @param object An NMRScaffold1D object.
#' @param getter Setter function.
#' @param slot.name Slot name as character.
#' @param value Value stored in specified slot.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name split_children
setGeneric("split_children", 
  function(object, setter, slot.name, value) standardGeneric("split_children")
)

#' @rdname split_children
#' @export
setMethod("split_children", "NMRScaffold1D", 
  function(object, setter, slot.name, value) {

  # If there are no children to split, then assign as normal
  if ( slot.name %in% slotNames(object) ) {
    slot(object, slot.name) <- value
    return(object)
  } 

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
#' Generic convenience method to access the id of any NMRScaffold1D object.
#' 
#' @param object An NMRScaffold1D object.
#' @param ... Additional arguments passed to inheriting methods.
#'
#' @name id
#' @export
setGeneric("id", 
  function(object, ...) standardGeneric("id")
)

#' @rdname id
#' @export
setMethod("id", "NMRScaffold1D", 
  function(object) object@id
)



#------------------------------------------------------------------------------
# Peaks

#---------------------------------------
#' Get object peaks
#' 
#' Generic convenience method to access the peak definitions of any
#' NMRScaffod1D object. If used on anything more complicated than an
#' NMRResonance1D, the function combines the peaks data frames of component
#' resonances.
#' 
#' @param object An NMRScaffold1D object.
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
setMethod("peaks", "NMRScaffold1D", 
  function(object, include.id = FALSE) {
  
  out <- combine_children(object, peaks, "peaks")
  if ( include.id & (nrow(out) > 0) ) {
    out <- cbind(id = object@id, out)
    colnames(out)[1] <- object@name
  }
  out

})

#---------------------------------------
#' Set object peaks
#' 
#' Generic convenience method to set the peak definitions of any NMRScaffold1D
#' object. This is primarily intended as an internal method, so use with
#' caution. Changing peak definitions after an object has been defined may have
#' unpredictable consequences.
#' 
#' @param object An NMRScaffold1D object.
#' @param value A data frame with a minimum of "position", "width", "height",
#'              and "fraction.gauss" columns. Peaks may be defined by one to
#'              three columns of "peak", "resonance", and "species" depending on
#'              the nature of the original object.
#' 
#' @name peaks-set
#' @export
setGeneric("peaks<-", 
  function(object, value) standardGeneric("peaks<-")
)

#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRScaffold1D",
  function(object, value) split_children(object, `peaks<-`, "peaks", value)
)


#------------------------------------------------------------------------------
# Couplings

#---------------------------------------
#' Get object couplings 
#' 
#' Generic convenience method to access the coupling definitions of any
#' NMRScaffold1D object. If used on anything more complicated than an
#' NMRResonance1D, the function combines the couplings data frames of component
#' resonances.
#' 
#' @param object An NMRScaffold1D object.
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
setMethod("couplings", "NMRScaffold1D", 
  function(object, include.id = FALSE) {

  out <- combine_children(object, couplings, "couplings")
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
#' Generic convenience method to access the bounds of any NMRScaffold1D object.
#' If used on anything more complicated than an NMRResonance1D, the function
#' combines the bounds data frames of component resonances. Whereas
#' lower_bounds() and upper_bounds() return a data.frame to match peaks(),
#' bounds() groups the lower and upper bounds into a list.
#' 
#' @param object An NMRScaffold1D object.
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
setMethod("bounds", "NMRScaffold1D", 
  function(object, include.id = FALSE) {
  
    list(lower = lower_bounds(object, include.id),
         upper = upper_bounds(object, include.id))

})

#' @rdname bounds
#' @export
setMethod("lower_bounds", "NMRScaffold1D", 
  function(object, include.id = FALSE) {

  out <- combine_children(object, lower_bounds, "lower.bounds")
  if ( include.id & (nrow(out) > 0) ) {
    out <- cbind(id = object@id, out)
    colnames(out)[1] <- object@name
  }
  out

})

#' @rdname bounds
#' @export
setMethod("upper_bounds", "NMRScaffold1D", 
  function(object, include.id = FALSE) {

  out <- combine_children(object, upper_bounds, "upper.bounds")
  if ( include.id & (nrow(out) > 0) ) {
    out <- cbind(id = object@id, out)
    colnames(out)[1] <- object@name
  }
  out

})



#==============================================================================>
# Convenience methods for bounds
#==============================================================================>



#---------------------------------------
#' Update bounds
#' 
#' This is an internal function used to update all the bounds of NMRScaffold1D object at once while ensuring that they are not widened if that is not desired. Do not use this function directly.
#' 
#' @param object An NMRScaffold1D object.
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
setMethod("update_bounds", "NMRScaffold1D", 
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
#' Set general bounds of an NMRScaffold1D object
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
#' @param object An NMRScaffold1D object.
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
setMethod("set_general_bounds", "NMRScaffold1D",
  function(object, position = NULL, height = NULL, width = NULL, 
           fraction.gauss = NULL, nmrdata = NULL, widen = FALSE) {
    
  # If object has children, just apply function recursively
  if ( "children" %in% slotNames(object) ) {

    object@children <- lapply(object@children, set_general_bounds, 
      position, height, width, fraction.gauss, nmrdata, widen
    )
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

    sfo1 <- get_parameter(nmrdata, 'sfo1', 'procs')
    width <- width * x.range * sfo1
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
  upper$fraction.gauss[upper$fraction.gauss > 0] <- 1

  # Ensuring that parameters are only widened if desired
  update_bounds(object, lower, upper, widen)
})


#------------------------------------------------------------------------------
#' Set offset bounds of an NMRScaffold1D object
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
#' @param object An NMRScaffold1D object.
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
setMethod("set_offset_bounds", "NMRScaffold1D",
  function(object, position = NULL, height = NULL, width = NULL, 
           relative = FALSE, widen = FALSE) {
    
  # If object has children, just apply function recursively
  if ( "children" %in% slotNames(object) ) {
    object@children <- lapply(object@children, set_offset_bounds, 
      position, height, width, relative, widen
    )
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
#' Set conservative bounds on an NMRScaffold1D object
#' 
#' A convenience function that sets reasonable bounds on an NMRScaffold1D
#' object. These bounds are assumed to be widely applicable to most simple NMR
#' data. Each set of bounds can be turned on or off as necessary. A slightly
#' better set of bounds can be selected if a reference NMRData1D object is
#' provided. If used on an anything more complicated than an NMRResonance1D,
#' the function propogates the bounds to component resonances.
#' 
#' @param object An NMRScaffold1D object.
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
setMethod("set_conservative_bounds", "NMRScaffold1D",
  function(object, position = TRUE,  height = TRUE, width = TRUE, 
           nmrdata = NULL, widen = FALSE) { 

  # First, do a single pass over general bounds with no reference
  if ( height )  gen.height <- c(0, Inf)
  else gen.height <- NULL

  if ( width ) gen.width <- c(0.003, 3)
  else gen.width <- NULL

  object <- set_general_bounds(object, height = gen.height, width = gen.width,
                               widen = widen)

  # Adding position offsets
  if ( position ) {
    object <- set_offset_bounds(object, position = c(-0.1, 0.1), widen = widen)
  }

  # If nmrdata is provided, add further constraints  
  if (! is.null(nmrdata) ) {
    
    if ( class(nmrdata) != 'NMRData1D' ) {
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
#' Set peak type of an NMRScaffold1D object
#' 
#' The peak type of an NMRScaffold1D object is governed by the fraction.gauss
#' parameter and can fluidly go from pure Lorentz to Gauss via the combined
#' Voigt lineshape. However, it can be computationally efficient to force a
#' specific peak type. This function provides a shortcut for doing so by
#' changing the fraction.gauss parameter as well as lower and upper bounds.
#' 
#' @param object An NMRScaffold1D object.
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
setMethod("set_peak_type", "NMRScaffold1D",
  function(object, peak.type) {
    
  # If object has children, just apply function recursively
  if ( "children" %in% slotNames(object) ) {
    object@children <- lapply(object@children, set_peak_type, peak.type) 
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
#  Initialization functions (generating parameter estimates based on data)
#==============================================================================>



#------------------------------------------------------------------------------
#' Initialize peak heights of an NMRScaffold1D object
#' 
#' Generates peak height estimates based on spectral data. If used on an
#' anything more complicated than an NMRResonance1D, the function propogates
#' the bounds to component resonances. Since estimates cannot be provided for
#' any peak outside the given data range, any such peaks are ignored. Whether
#' or not warning/error messages are generated when that occurs is specified by
#' the exclusion.notification parameter. Similarly, exclusion.level provided
#' options to omit the whole resonance/species if a peak is found to be outside
#' the data bounds.
#' 
#' At this point, there is just one approach: take peak height as the intensity
#' of the data at the current position of the peak. There are plans to develop
#' more sophisticated approaches in the future.
#' 
#' @param object An NMRScaffold1D object.
#' @param nmrdata An NMRData1D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A new object with modified peak heights.
#' 
#' @name initialize_heights
#' @export
setGeneric("initialize_heights", 
  function(object, ...) standardGeneric("initialize_heights"))

#' @rdname initialize_heights
#' @export
setMethod("initialize_heights", "NMRScaffold1D",
  function(object, nmrdata) {

  # Checking nmrdata
  if ( class(nmrdata) != 'NMRData1D' ) {
    err <- '"nmrdata" must be a valid NMRData1D object.'
    stop(err)
  }

  # Checking that data captures all defined peaks
  d <- processed(nmrdata)
  peaks <- peaks(object) 

  logic <- (peaks$position < min(d$direct.shift)) | 
           (peaks$position > max(d$direct.shift))

  if ( any(logic) ) {
    err <- "nmrdata must span all peak positions (%s are out of bounds)"
    err <- sprintf(err, paste(peaks$id[logic], collapse = ", "))
    stop(err)
  }

  # Building an interpolating function betwewn chemical shift and intensity
  f <- approxfun(d$direct.shift, Re(d$intensity))

  # Generating heights from interpolation
  peaks$height <- f(peaks$position)

  # Updating
  peaks(object) <- peaks
  
  object
})



#========================================================================>
#  Lineshape and area calculations
#========================================================================>



#------------------------------------------------------------------------
#' Generate lineshape function
#' 
#' This is primarily an internal method that outputs a function (or a tbl_df
#' data frame of functions), where each function outputs spectral intensity
#' data given a vector input of chemical shifts.
#' 
#' @param object An NMRScaffold1D object.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmroptions$direct$sf = ..., but an override can be
#'           provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single function, FALSE to output a data frame of functions
#'                  that correspond to individual peaks.
#' @param include.id TRUE to include id outer column if outputting data frame.
#' @param components 'r/i' to output both real and imaginary data, 'r' to output
#'                   only real and 'i' to output only imaginary.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A function or tbl_df data frame of functions where each function
#'         outputs spectral intensity data given a vector input of chemical
#'         shifts. In the latter case, the functions are stored in a list column
#'         called f.
#' 
#' @name f_lineshape
#' @export
setGeneric("f_lineshape", 
  function(object, ...) standardGeneric("f_lineshape")
)

#' @rdname f_lineshape
#' @export
setMethod("f_lineshape", "NMRScaffold1D",
  function(object, sf = nmroptions$direct$sf, sum.peaks = TRUE, 
           include.id = FALSE, components = 'r/i') {

    # Checking to make sure that sweep frequency is defined
    err <- '"sf" must be provided as input or set using nmroptions$direct$sf'
    if ( is.null(sf) ) stop(err)

    # Defining which components to return
    return.r <- grepl('r', tolower(components))
    return.i <- grepl('i', tolower(components))

    err <- '"components" must have at least one of either "r" or "i"'
    if ( return.r && return.i ) f_out <- function(y) {y}
    else if ( return.r ) f_out <- function(y) {Re(y)}
    else if ( return.i ) f_out <- function(y) {Im(y)}
    else stop(err)

    columns <- c('position', 'width', 'height', 'fraction.gauss')
    peaks <- peaks(object)
    parameters <- as.matrix(peaks[, columns])

    # Converting peak width to ppm
    parameters[, 2] <- parameters[, 2]/sf

    # If peaks are to be summed, just feed all parameters into the Rcpp function
    if ( sum.peaks ) {
      out <- function(x) {
        
        n <- as.integer(length(x))
        y <- .Call("eval_1d_wrapper",        
          x = as.double(x),
          y = as.double(rep(0, n*2)),
          knots = as.double(0),
          p = as.double(as.vector(t(parameters))),
          n = n,
          nl = as.integer(length(parameters)),
          nb = as.integer(0),
          np = as.integer(0),
          nk = as.integer(0)
        )

        f_out(cmplx1(r = y[1:n], i = y[(n+1):(n*2)]))
      }
    } 
    # Otherwise, generate a tbl_df data frame
    else {
      out <- as_tibble(peaks[, which(! colnames(peaks) %in% columns)])
      if ( include.id && (nrow(out) > 0) ) {
        if ( 'resonance' %in% colnames(out) ) {
          out <- cbind(species = object@id, out)
        }
        else {
          out <- cbind(resonance = object@id, out)
        }
      }

      parameters <- split(parameters, 1:nrow(parameters))
      
      # Generating a list of functions, each with their parameters enclosed
      functions <- lapply(parameters, function (p) {
        function(x) {
          n <- as.integer(length(x))
          y <- .Call("eval_1d_wrapper",        
            x = as.double(x),
            y = as.double(rep(0, n*2)),
            knots = as.double(0),
            p = as.double(as.vector(p)),
            n = n,
            nl = as.integer(length(p)),
            nb = as.integer(0),
            np = as.integer(0),
            nk = as.integer(0)
          )

          f_out(cmplx1(r = y[1:n], i = y[(n+1):(n*2)]))
        }
      })

      # Adding functions as a column
      out$f <- functions
    }

    out
  })



#------------------------------------------------------------------------
#' Calculate peak lineshape values
#' 
#' Calculated peak intensity values over a set of chemical shifts.
#' 
#' @param object An NMRScaffold1D object.
#' @param direct.shift Vector of chemical shift data in ppm.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmroptions$direct$sf = ..., but an override can be
#'           provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single set of values, FALSE to output a data frame of values
#'                  that correspond to individual peaks.
#' @param sum.baseline TRUE to add baseline to every peak, if one is defined.
#'                     FALSE to exclude baseline. 
#' @param include.id TRUE to include id as outer column if outputting data
#'                   frame.
#' @param components 'r/i' to output both real and imaginary data, 'r' to output
#'                   only real and 'i' to output only imaginary.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A vector of spectral intensity data or a tibble with columns
#'         "resonance" (optional), "peak", "direct.shift", and "intensity".
#' 
#' @name values
#' @export
setGeneric("values", 
  function(object, ...) standardGeneric("values")
)

#' @rdname values
#' @export
setMethod("values", "NMRScaffold1D",
  function(object, direct.shift, sf = nmroptions$direct$sf, sum.peaks = TRUE, 
           sum.baseline = FALSE, include.id = FALSE, components = 'r/i') {

  # Generating baseline if necessary
  if ( sum.baseline && (class(object) == 'NMRFit1D') ) {
    f <- f_baseline(object, components)
    baseline <- f(direct.shift)
  } else {
    baseline <- rep(0, length(direct.shift))
  }

  # Output depends on whether peaks are summed or not
  if ( sum.peaks ) {
    # Get function
    f <- f_lineshape(object, sf, sum.peaks, components)

    # And apply it to specified chemical shifts
    f(direct.shift) + baseline
  } 
  else {
    # Get data frame of functions
    d <- f_lineshape(object, sf, sum.peaks, include.id, components)

    # Defining function that generates necessary data frame
    f <- function(d) {
      intensity <- d$f[[1]](direct.shift) + baseline
      data.frame(direct.shift = direct.shift, 
                 intensity = vec_cast(intensity, complex())) 
    }

    # And apply it for every peak
    descriptors <- colnames(d)[1:(ncol(d)-1)]
    d %>% 
      group_by_at(descriptors) %>% 
      do( f(.) ) %>%
      ungroup()
  }
})



#------------------------------------------------------------------------
#' Calculate peak areas
#' 
#' Calculate total peak areas based on peak parameters.
#' 
#' @param object An NMRScaffold1D object.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmroptions$direct$sf, but an override can be
#'           provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single area, FALSE to output a data frame of peak area
#'                  values.
#' @param include.id TRUE to include id as outer column if outputting data
#'                   frame.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A single overall area or a data frame of areas with columns
#'         "resonance" (optional), "peak", and "area".
#' 
#' @name areas
#' @export
setGeneric("areas", 
  function(object, ...) standardGeneric("areas")
)

#' @rdname areas 
#' @export
setMethod("areas", "NMRScaffold1D",
  function(object, sf = nmroptions$direct$sf, sum.peaks = TRUE, 
           include.id = FALSE, components = 'r/i') {

  # Defining area function
  f <- function(position, width, height, fraction.gauss) {
    # If fraction is 0, treat as Lorentz
    if ( fraction.gauss == 0 ) {
      pi*width*height
    }
    # If fraction is 1, treat as Gauss
    else if ( fraction.gauss == 1) {
      sqrt(2*pi)*width*height
    }
    # Else, proceed as Voigt
    else {
      l.width <- width
      g.width <- width*fraction.gauss/(1 - fraction.gauss)
      Re(sqrt(2*pi)*g.width*height /
         Faddeeva_w(complex(im = l.width)/(sqrt(2)*g.width)))
    }
  }

  # Calculating areas
  areas <- peaks(object, include.id)

  areas <- areas %>%
    group_by_all() %>%
    summarize(area = f(position, width, height, fraction.gauss)) %>%
    ungroup() %>%
    select(-position, -width, -height, -fraction.gauss) %>%
    as.data.frame()

  # Sum if necessary
  if ( sum.peaks ) sum(areas$area)
  else areas
})
