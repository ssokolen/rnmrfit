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
#' @slot dimensions A list of NMRResonance1D object.
#' 
#' @name NMRResonance2D-class
#' @export
NMRResonance2D <- setClass("NMRResonance2D",
  contains = "NMRScaffold2D",
  slots = c(
    name = "character",
    id = "character",
    dimensions = "list"
  ),
  prototype = prototype(
    name = 'resonance',
    id = 'resonance'
  )
)



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

  # If an id is not provided, generate one
  if ( is.null(id) ) {

    if ( direct@id != indirect@id ) {
      id <- paste(direct@id, indirect@id, sep = ' / ')
    } else {
      id <- direct@id
  }
}

  #---------------------------------------
  # Generate object

  dimensions <- list(direct = direct, indirect = indirect)
  new('NMRResonance2D', id = id, dimensions = dimensions)
}



#==============================================================================>
#  Bounds
#==============================================================================>



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


