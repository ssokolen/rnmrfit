# Definition of a class structure for 2D resonance data.



#==============================================================================>
#  NMRResonance2D -- peak description 
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR resonance.
#' 
#' Essentially, this class is used to define coupling relationships to group
#' individual peaks into resonances. A resonance is a series of peaks that are
#' constrained by a set of relations between their position, width, and,
#' fraction.gauss parameters as well as overall areas. A 2D resonance is simply
#' a collection of 1D resonances in the direct and indirect dimension with
#' separate set of coupling constraints.
#' 
#' @slot dimensions An NMRResonance1D object.
#' @slot indirect An NMRResonance1D object.
#' 
#' @name NMRResonance2D-class
#' @export
NMRResonance2D <- setClass("NMRResonance2D",
  contains = 'NMRScaffold2D',
  slots = c(
    direct = 'NMRResonance1D',
    indirect = 'NMRResonance1D'
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRResonance2D validity test
#'
validNMRResonance2D <- function(object) {

  direct <- direct(object)
  indirect <- indirect(object)

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Both direct and indirect objects must be NMRResonance1D objects
  logic1 <- class(direct) == 'NMRResonance1D'
  logic2 <- class(indirect) == 'NMRResonance1D'

  if (! (logic1 && logic2) ) {

    valid <- FALSE
    new.err <- paste('"direct" and "indirect" components must be valid',
                     'NMRResonance1D objects.')
    err <- c(err, new.err)
  }   
  
  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRResonance2D", validNMRResonance2D)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRResonance2D object based on simplified peak list
#' 
#' Generates an NMRResonance2D object by combining direct and indirect peak
#' definitions. These peak definitions can be provided either with a string
#' that will be parsed by \code{parse_peaks_1d()} or a previously generated
#' NMRResonance1D object.
#' 
#' @param direct.peaks An NMRResonance1D object, numeric vector of singlet
#'                     chemical shifts, or a character string specifying
#'                     multiplets of the form "3 d 1.2". See ?parse_peaks_1d for
#'                     more information.
#' @param indirect.peaks An NMRResonance1D object, numeric vector of singlet
#'                       chemical shifts, or a character string specifying
#'                       multiplets of the form "3 d 1.2". See ?parse_peaks_1d
#'                       for more information.
#' @param direct.sf Sweep frequency (MHz) in the direct dimension -- needed to
#'                  convert coupling constants from Hz to ppm. In most cases, it
#'                  is recommended to set a single default value using
#'                  nmroptions$direct$sf  = ..., but an override can be
#'                  provided here.
#' @param indirect.sf Sweep frequency (MHz) in the indirect dimension -- needed
#'                    to convert coupling constants from Hz to ppm. In most
#'                    cases, it is recommended to set a single default value
#'                    using nmroptions$indirect$sf  = ..., but an override
#'                    can be provided here.
#' @param id A string specifying resonance name. If left empty, a name is
#'           automatically generated from the peaks argument.
#' @param width Initial estimate of peak width (in Hz). For Voigt lineshapes,
#'              this value is taken as the Lorentzian component, with the
#'              Gaussian component calculated from
#'              peak.width*frac.guass/(1-frac.gauss).
#' @param fraction.gauss Fraction of overall peak width that corresponds to a
#'                       Gaussian lineshape. A value of 0 corresponds to a
#'                       Lorentz peak whereas a value of 1 corresponds to a
#'                       Gaussian peak. Values in between 0 and 1 are modelled
#'                       as a Voigt lineshape but the specific value of
#'                       frac.gauss does not have a physical interpretation.
#' @param position.leeway A fraction specifying how tightly enforced the
#'                        coupling constraints on peak positions, should be.
#'                        E.g. coupling.leeway = 0 specifies that the j coupling
#'                        constant is exact, whereas couping.leeway = 0.1
#'                        specifies that the coupling constant may differ by +/-
#'                        10 percent.
#' @param width.leeway Similar to position.leeway but for peak widths.
#'                     Determines how strictly equal peak widths for all coupled
#'                     peaks are enforced.
#' @param area.leeway Similar to position.leeway but for peak areas. Determines
#'                    how strictly the coupling area ratios are enforced.
#' 
#' @return An NMRResonance2D object.
#' 
#' @export
nmrresonance_2d <- function(direct.peaks, indirect.peaks, 
                            direct.sf = nmroptions$direct$sf, 
                            indirect.sf = nmroptions$indirect$sf, 
                            id = NULL, width = 1, fraction.gauss = 0, 
                            position.leeway = 0, area.leeway = 0, 
                            width.leeway = 0) {

  #---------------------------------------
  # First, generating the direct/indirect components

  # If the direct.peaks component is already an NMRResonance1D object,
  # then there is nothing to do, otherwise, pass input into nmrresonance_1d
  if ( 'NMRResonance1D' %in% class(direct.peaks) ) {
    direct <- direct.peaks
  } else {
    direct <- nmrresonance_1d(direct.peaks, sf = direct.sf, id = id,
                              width = width, fraction.gauss = fraction.gauss,
                              position.leeway = position.leeway,
                              area.leeway = area.leeway, 
                              width.leeway = width.leeway)
  }

  # Same for indirect
  if ( 'NMRResonance1D' %in% class(indirect.peaks) ) {
    indirect <- indirect.peaks
  } else {
    indirect <- nmrresonance_1d(indirect.peaks, sf = indirect.sf, id = id,
                                width = width, fraction.gauss = fraction.gauss,
                                position.leeway = position.leeway,
                                area.leeway = area.leeway, 
                                width.leeway = width.leeway)
  }

  #---------------------------------------
  # Aligning the ids

  # If an id is provided apply it directly to both components
  if (! is.null(id) ) {
    direct@id <- id
    indirect@id <- id
  }
  # If the id is not provided, and individual ids are different, combine them.
  else if ( direct@id != indirect@id ) {
    id <- paste(direct@id, indirect@id, sep = ' / ')
    direct@id <- id
    indirect@id <- id
}

  #---------------------------------------
  # Generate object

  new('NMRResonance2D', direct = direct, indirect = indirect)
}



#==============================================================================>
#  Helper functions
#==============================================================================>



#---------------------------------------
#' Combine direct and indirect dimensions
#' 
#' This is an internal function used for all getter functions that output a
#' data.frame object. Essentially, the getter is passed on to the direct and
#' indirect components, a dimension column is added and the resulting objects
#' are stitched together.
#' 
#' @param object An NMRResonance2D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name combine_dimensions
setGeneric("combine_dimensions", 
  function(object, ...) standardGeneric("combine_dimensions")
)

#' @rdname combine_dimensions
setMethod("combine_dimensions", "NMRResonance2D", 
  function(object, getter, ...) {
    direct <- data.frame(dimension = "direct", 
                         getter(object@direct, ...))
    indirect <- data.frame(dimension = "indirect", 
                           getter(object@indirect, ...))
    rbind(direct, indirect)
})

#---------------------------------------
#' Split direct and indirect dimensions
#' 
#' This is an internal function used for all setter functions that input a
#' data.frame object. Essentially, the input value is first split based on a
#' "dimension" column, with the setter being passed on to the direct and
#' indirect components.
#' 
#' @param object An NMRResonance2D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name split_dimensions
setGeneric("split_dimensions", 
  function(object, setter, value) standardGeneric("split_dimensions")
)

#' @rdname split_dimensions
#' @export
setMethod("split_dimensions", "NMRResonance2D", 
  function(object, setter, value) {

    # First the input must be a data.frame of some sort
    err <- 'Input value must be a data.frame type object.'
    if (! 'data.frame' %in% class(value) ) stop(err)

    # Second the input must have "dimension" column
    err <- 'Input data.frame must have a "dimension" column.'
    if (! 'dimension' %in% colnames(value) ) stop(err)

    # Third, the dimension column must only contain direct and indirect values
    err <- 'The "dimension" column must only contain "direct" or "indirect".'
    entries <- sort(unique(as.character(value$dimension)))
    if (! identical(entries, c('direct', 'indirect')) ) stop(err)

    # If all of the above is met, then split components
    direct <- filter(value, dimension == 'direct') %>% select(-dimension)
    object@direct <- setter(object@direct, direct)

    indirect <- filter(value, dimension == 'indirect') %>% select(-dimension)
    object@indirect <- setter(object@indirect, indirect)
    
    validObject(object)
    object
})



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>



#' @rdname id
#' @export
setMethod("id", "NMRResonance2D", 
  function(object) object@direct@id
)



#' @rdname id-set
#' @export
setReplaceMethod("id", "NMRResonance2D",
  function(object, value) {
    id <- as.character(value)
    object@direct@id <- id
    object@indirect@id <- id
    validObject(object)
    object 
  })



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
# Direct and indirect components



#' @rdname direct
#' @export
setMethod("direct", "NMRResonance2D", 
  function(object) object@direct
)



#' @rdname indirect
#' @export
setMethod("indirect", "NMRResonance2D", 
  function(object) object@indirect
)



#------------------------------------------------------------------------------
# Id



#' @rdname id
#' @export
setMethod("id", "NMRResonance2D", 
  function(object) object@direct@id
)



#' @rdname id-set
#' @export
setReplaceMethod("id", "NMRResonance2D",
  function(object, value) {
    id <- as.character(value)
    object@direct@id <- id
    object@indirect@id <- id
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Peaks



#' @rdname peaks
#' @export
setMethod("peaks", "NMRResonance2D", 
  function(object, include.id = FALSE) {
    combine_dimensions(object, peaks, include.id = include.id)
})



#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRResonance2D",
  function(object, value) {
    split_dimensions(object, `peaks<-`, value)
})



#' @rdname update_peaks
setMethod("update_peaks", "NMRResonance2D",
  function(object, peaks, exclusion.level = nmroptions$exclusion$level,
           exclusion.notification = nmroptions$exclusion$notification) {
    split_dimensions(object, update_peaks, peaks, 
                     exclusion.level, exclusion.notification)
})



#==============================================================================>
#  Bounds
#==============================================================================>



#' @rdname set_general_bounds
#' @export
setMethod("set_general_bounds", "NMRResonance2D",
  function(object, dimension = "both", position = NULL, height = NULL, 
           width = NULL, fraction.gauss = NULL, nmrdata = NULL, widen = FALSE) {

  # Checking dimension argument
  err <- 'dimension must be one of "direct", "indirect", or "both"'
  if (! dimension %in% c("direct", "indirect", "both") ) stop(err)

  # Checking data
  if (! is.null(nmrdata) ) {

    if ( class(nmrdata) != 'NMRData2D' ) {
      err <- '"nmrdata" must be a valid NMRData2D object.'
      stop(err)
    }
    else {
      validObject(nmrdata)
    }
  }

  # Adding direct dimension bounds, if necessary
  if ( dimension %in% c('direct', 'both') ) {

    if (! is.null(nmrdata) ) {
      columns <- c('direct.shift', 'intensity')
      d <- new("NMRData1D", processed = processed(nmrdata)[, columns])
    } else {
      d <- NULL
    }

    object@direct <- set_general_bounds(object@direct, position, height,
                                        width, fraction.gauss, d, widen)
  }

  # And then the indirect
  if ( dimension %in% c('indirect', 'both') ) {

    if (! is.null(nmrdata) ) {
      columns <- c('indirect.shift', 'intensity')
      d <- new("NMRData1D", processed = processed(nmrdata)[, columns])
    } else {
      d <- NULL
    }

    object@indirect <- set_general_bounds(object@indirect, position, height,
                                          width, fraction.gauss, d, widen)
  }

  validObject(object)
  object
})



#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRResonance2D",
  function(object, dimension = "both", position = NULL, height = NULL, 
           width = NULL, relative = FALSE, widen = FALSE) {

  # Checking dimension argument
  err <- 'dimension must be one of "direct", "indirect", or "both"'
  if (! dimension %in% c("direct", "indirect", "both") ) stop(err)

  # Adding direct dimension bounds, if necessary
  if ( dimension %in% c('direct', 'both') ) {
    object@direct <- set_offset_bounds(object@direct, position, height,
                                       width, relative, widen)
  }

  # And then the indirect
  if ( dimension %in% c('indirect', 'both') ) {
    object@indirect <- set_offset_bounds(object@indirect, position, height,
                                         width, relative, widen)
  }

  validobject(object)
  object
})



#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRResonance2D",
  function(object, position = TRUE,  height = TRUE, width = TRUE, 
           nmrdata = NULL, widen = FALSE) { 

  # Checking dimension argument
  err <- 'dimension must be one of "direct", "indirect", or "both"'
  if (! dimension %in% c("direct", "indirect", "both") ) stop(err)

  # Checking data
  if (! is.null(nmrdata) ) {

    if ( class(nmrdata) != 'NMRData2D' ) {
      err <- '"nmrdata" must be a valid NMRData2D object.'
      stop(err)
    }
    else {
      validObject(nmrdata)
    }
  }

  # Adding direct dimension bounds, if necessary
  if ( dimension %in% c('direct', 'both') ) {

    if (! is.null(nmrdata) ) {
      columns <- c('direct.shift', 'intensity')
      d <- new("NMRData1D", processed = processed(nmrdata)[, columns])
    } else {
      d <- NULL
    }

    object@direct <- set_conservative_bounds(
      object@direct, position, height, width, fraction.gauss, d, widen
    )
  }

  # And then the indirect
  if ( dimension %in% c('indirect', 'both') ) {

    if (! is.null(nmrdata) ) {
      columns <- c('indirect.shift', 'intensity')
      d <- new("NMRData1D", processed = processed(nmrdata)[, columns])
    } else {
      d <- NULL
    }

    object@indirect <- set_conservative_bounds(
      object@indirect, position, height, width, fraction.gauss, d, widen
    )
  }

  validObject(object)
  object
})


#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRResonance2D",
  function(object, peak.type) {

  object@direct <- set_peak_type(object@direct, peak.type)
  object@indirect <- set_peak_type(object@indirect, peak.type)

  validObject(object)
  object
})


