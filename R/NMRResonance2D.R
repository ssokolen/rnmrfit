# Definition of a class structure for 2D resonance data.



#==============================================================================>
#  NMRResonance2D -- peak description 
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR resonance.
#' 
#' Essentially, this class is used to define coupling relationships to group
#' individual peaks into resonances. A resonances is a series of peaks that are
#' constrained by a set of relations between their position, width, and,
#' fraction.gauss parameters as well as overall areas. A 2D resonance is simply
#' a collection of 1D resonances in the direct and indirect dimension with
#' separate set of coupling constraints.
#' 
#' @slot direct An NMRResonance1D object.
#' @slot indirect An NMRResonance1D object.
#' 
#' @name NMRResonance2D-class
#' @export
NMRResonance2D <- setClass("NMRResonance2D",
  contains = 'NMRScaffold2D',
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRResonance2D validity test
#'
validNMRResonance2D <- function(object) {

  direct <- object@direct
  indirect <- object@indirect

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

  # If an id is provided apply it directly to both component
  }

  #---------------------------------------
  # Generate object

  new('NMRResonance2D', direct = direct, indirect = indirect)
}
