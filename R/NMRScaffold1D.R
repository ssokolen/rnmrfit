# Definition of a super-class for 1D resonance data.

#' @import Rcpp
#' @useDynLib rnmrfit
NULL


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
#' Display any NMRScaffold1D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRScaffold1D", 
  function(object) {

    # Generating compiled data frames
    id <- ifelse( 'id' %in% slotNames(object), sprintf('(%s)', id(object)), '')
    peaks <- peaks(object)
    bounds <- bounds(object)
    couplings <- couplings(object)

    # Id (if it exists)
    cat(sprintf('An object of %s class %s\n\n', class(object), id))

    # Peaks
    cat('Peaks:\n\n')
    print(peaks)
    cat('\n')

    # Bounds
    columns <- c('position', 'width', 'height', 'fraction.gauss')

    lower <- unlist(bounds$lower[ , columns])
    upper <- unlist(bounds$upper[ , columns])
    
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

#---------------------------------------
#' Set object id
#' 
#' Generic convenience method to set the id of any NMRScaffold1D object.
#' 
#' @param object An NMRScaffold1D object.
#' @param value New id (converted to character).
#'
#' @name id-set
#' @export
setGeneric("id<-", 
  function(object, value) standardGeneric("id<-")
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


#------------------------------------------------------------------------------
#' Update peak parameters
#' 
#' This is an internal function whose role is to update existing peak
#' parameters, while accounting for exclusion criteria and generating relevant
#' errors/messages. Any peak not updated is excluded.
#' 
#' @param object An NMRScaffold1D object.
#' @param peaks A data frame with "position", "width", "height", and
#'              "fraction.gauss" columns. Peaks may be defined by one to three
#'              columns of "peak", "resonance", and "species" depending on the
#'              nature of the original object.
#' @param exclusion.level A string specifying what to do when peaks are found to
#'                        fall outside of the data range: either 'species' to
#'                        exclude the whole species to which the offending peak
#'                        belongs, 'resonance' to exclude the resonance to which
#'                        the offending peak belongs, or 'peak' to exclude just
#'                        the peak itself.
#' @param exclusion.notification A function specifying how to report when peaks
#'                               are found to be outside the data range: 'none'
#'                               to ignore, 'message' to issue a message,
#'                               'warning' to issue a warning, and 'stop' to
#'                               issue an error.
#' @inheritParams methodEllipse
#' 
#' @return A new object with modified peak parameters.
#' 
#' @name update_peaks
setGeneric("update_peaks", 
  function(object, ...) standardGeneric("update_peaks")
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

#---------------------------------------
#' Set object couplings
#' 
#' Generic convenience method to set the coupling definitions of an
#' NMRResonance1D object.
#' 
#' @param object An NMRResonance1D object.
#' @param value A data frame with "id.1", "id.2", "position.difference", and
#'              "area.ratio" columns.
#' 
#' @name couplings-set
#' @export
setGeneric("couplings<-", 
  function(object, value) standardGeneric("couplings<-")
)



#------------------------------------------------------------------------------
# Bounds

#---------------------------------------
#' Get object bounds
#' 
#' Generic convenience method to access the bounds of any NMRScaffold1D object.
#' If used on an anything more complicated than an NMRResonance1D, the function
#' combines the bounds data frames of component resonances.
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

#---------------------------------------
#' Set object bounds
#' 
#' Generic convenience method to set the bounds of an NMRResonance1D object.
#' 
#' @param object An NMRResonance1D object.
#' @param value A list with "lower" and "upper" elements, each containing a data
#'              frame with "peak", "position", "width", "height", and
#'              "fraction.gauss" columns.
#' 
#' @name bounds-set
#' @export
setGeneric("bounds<-", 
  function(object, value) standardGeneric("bounds<-")
)



#==============================================================================>
# Convenience methods for bounds
#==============================================================================>



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
#' @param exclusion.level A string specifying what to do when peaks are found to
#'                        fall outside of the data range: either 'species' to
#'                        exclude the whole species to which the offending peak
#'                        belongs, 'resonance' to exclude the resonance to which
#'                        the offending peak belongs, or 'peak' to exclude just
#'                        the peak itself.
#' @param exclusion.notification A function specifying how to report when peaks
#'                               are found to be outside the data range: 'none'
#'                               to ignore, 'message' to issue a message,
#'                               'warning' to issue a warning, and 'stop' to
#'                               issue an error.
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
  function(object, nmrdata, exclusion.level = nmroptions$exclusion$level,
           exclusion.notification = nmroptions$exclusion$notification) {

  # Checking nmrdata
  if ( class(nmrdata) != 'NMRData1D' ) {
    err <- '"nmrdata" must be a valid NMRData1D object.'
    stop(err)
  }
  else {
    validObject(nmrdata)
  }

  # Building an interpolating function betwewn chemical shift and intensity
  d <- processed(nmrdata)
  f <- approxfun(d$direct.shift, Re(d$intensity))

  # Excluding peaks that are outside the data frame
  peaks <- peaks(object) 
  logic <- (peaks$position > min(d$direct.shift)) & 
           (peaks$position < max(d$direct.shift))

  peaks <- peaks[logic, ]

  # Generating heights from interpolation
  peaks$height <- f(peaks$position)

  # Updating
  object.2 <- update_peaks(object, peaks, exclusion.level = exclusion.level,
                           exclusion.notification = exclusion.notification)
  peaks.2 <- peaks(object.2)
  
  all.columns <- colnames(peaks)
  data.columns <- c('position', 'width', 'height', 'fraction.gauss')
  id.columns <- all.columns[! all.columns %in% data.columns]

  peaks <- peaks(object)
  ids <- apply(peaks[, id.columns], 1, paste, collapse = '-')
  ids.2 <- apply(peaks.2[, id.columns], 1, paste, collapse = '-')

  peaks[ids %in% ids.2, ] <- peaks.2
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
#' @param object An NMRScaffold1D or NMRScaffold2d object.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. A single value can be used for 1D and homonuclear 2D
#'           experiments or two values for the direct and indirect dimensions,
#'           in that order. this can be  In most cases, it is recommended to set
#'           a single default value using, for example, nmroptions$direct$sf =
#'           ..., but an override can be provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single function, FALSE to output a data frame of functions
#'                  that correspond to individual peaks.
#' @param include.id TRUE to include id outer column if outputting data frame.
#' @param components A string specifying the real/imaginary components to
#'                   output, e.g. 'r/i' or 'rr/ri/ir/ii'. Use the symbols 'r',
#'                   'i' for 1D experiments and 'rr', 'ri', 'ir', 'ii' for 2D
#'                   experiments. The symbols can be separated by spaces or
#'                   practically any other characters, so 'r i' is interpreted
#'                   the same as 'r/i'. The total number of components will
#'                   determine the final data type of the output, with 1
#'                   component corresponding to real values, 2 components to
#'                   cmplx1d objects (similar to the normal complex type), and 3
#'                   or 4 components to cmplx2d hypercomplex objects.
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
    components <- rev(sort(strsplit(components, '[^ri]+', perl = TRUE)[[1]]))

    # Checking 1D components
    err <- paste('"component" argument must consist of one-character codes',
                 'and possibly a separator, e.g., "r/i" or "r"')
    if ( any(! components %in% c('r', 'i')) ) stop(err)

    return.r <- "r" %in% components
    return.i <- "i" %in% components

    if ( return.r && return.i ) f_out <- function(y) {y}
    else if ( return.r ) f_out <- function(y) {Re(y)}
    else if ( return.i ) f_out <- function(y) {Im(y)}

    columns <- c('position', 'width', 'height', 'fraction.gauss')
    peaks <- peaks(object, include.id = include.id)
    parameters <- as.matrix(peaks[, columns])

    # Converting peak width to ppm
    parameters[, 2] <- parameters[, 2]/sf

    # Defining function generator based on arbitrary subset of parameters
    f_gen <- function(p) {
      force(p)
      function(x) {
        y <- matrix(0, nrow = length(x), ncol = 2)
        .Call('_rnmrfit_lineshape_1d', PACKAGE = 'rnmrfit', x, y, p)
        f_out(cmplx1(r = y[,1], i = y[,2]))
      }
    }


    # If peaks are to be summed, just feed all parameters into the Rcpp function
    if ( sum.peaks ) {
      p <- as.vector(t(parameters))
      out <- f_gen(p)
    } 
    # Otherwise, generate a tbl_df data frame
    else {
      out <- as_tibble(peaks[, which(! colnames(peaks) %in% columns)])
      parameters <- split(parameters, 1:nrow(parameters))
      
      # Generating a list of functions, each with their parameters enclosed
      functions <- lapply(parameters, function (p) {
        p <- as.vector(t(p))
        f_gen(p)
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
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @param direct.shift Vector of chemical shift data in ppm.
#' @param indirect.shift Vector of chemical shift data in ppm.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. A single value can be used for 1D and homonuclear 2D
#'           experiments or two values for the direct and indirect dimensions,
#'           in that order. this can be  In most cases, it is recommended to set
#'           a single default value using, for example, nmroptions$direct$sf =
#'           ..., but an override can be provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single set of values, FALSE to output a data frame of values
#'                  that correspond to individual peaks.
#' @param sum.baseline TRUE to add baseline to every peak, if one is defined.
#'                     FALSE to exclude baseline.
#' @param include.id TRUE to include id as outer column if outputting data
#'                   frame.
#' @param components A string specifying the real/imaginary components to
#'                   output, e.g. 'r/i' or 'rr/ri/ir/ii'. Use the symbols 'r',
#'                   'i' for 1D experiments and 'rr', 'ri', 'ir', 'ii' for 2D
#'                   experiments. The symbols can be separated by spaces or
#'                   practically any other characters, so 'r i' is interpreted
#'                   the same as 'r/i'. The total number of components will
#'                   determine the final data type of the output, with 1
#'                   component corresponding to real values, 2 components to
#'                   cmplx1d objects (similar to the normal complex type), and 3
#'                   or 4 components to cmplx2d hypercomplex objects.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A vector of spectral intensity data or a tibble with peak
#'         identifiers.
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
    f <- f_lineshape(object, sf, sum.peaks, include.id, components)

    # And apply it to specified chemical shifts
    f(direct.shift) + baseline
  } 
  else {
    # Get data frame of functions
    d <- f_lineshape(object, sf, sum.peaks, include.id, components)

    # Defining function that generates necessary data frame
    f <- function(g) {
      tibble(direct.shift = direct.shift, 
             intensity = (g[[1]](direct.shift)) + baseline) %>%
      unpack('intensity')
    }

    # Note that the unpack/pack functions are used to avoid bind_row errors
    
    # And apply them to every peak
    d %>%
      group_by_if(function(x) {!is.list(x)}) %>% 
      do(f(.$f) ) %>% 
      pack('intensity')
  }
})



#------------------------------------------------------------------------
#' Calculate peak areas
#' 
#' Calculate total peak areas based on peak parameters.
#' 
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single area, FALSE to output a data frame of peak area
#'                  values.
#' @param include.id TRUE to include id as outer column if outputting data
#'                   frame.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A single overall area or a data frame of areas with peak identifier
#'         and area columns.
#' 
#' @name areas
#' @export
setGeneric("areas", 
  function(object, ...) standardGeneric("areas")
)

#' @rdname areas 
#' @export
setMethod("areas", "NMRScaffold1D",
  function(object, sum.peaks = TRUE, include.id = FALSE) {

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
