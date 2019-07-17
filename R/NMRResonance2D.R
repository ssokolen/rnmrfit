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
#' @param sf Sweep frequency (MHz) in the direct dimension -- needed to convert
#'           coupling constants from Hz to ppm. A single value can be used for
#'           homonuclear experiments or two values for the direct and indirect
#'           dimensions, in that order. In most cases, it is recommended to set
#'           a single default value using, for example, nmroptions$direct$sf  =
#'           ..., but an override can be provided here.
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
                            sf = c(nmroptions$direct$sf, 
                                   nmroptions$indirect$sf), 
                            id = NULL, width = 1, fraction.gauss = 0, 
                            position.leeway = 0, area.leeway = 0, 
                            width.leeway = 0) {

  if ( length(sf) == 1 ) sf <- c(sf, sf)
 
  #---------------------------------------
  # First, generating the direct/indirect components

  # If the direct.peaks component is already an NMRResonance1D object,
  # then there is nothing to do, otherwise, pass input into nmrresonance_1d
  if ( 'NMRResonance1D' %in% class(direct.peaks) ) {
    direct <- direct.peaks
  } else {
    direct <- nmrresonance_1d(direct.peaks, sf = sf[1], id = id,
                              width = width, fraction.gauss = fraction.gauss,
                              position.leeway = position.leeway,
                              area.leeway = area.leeway, 
                              width.leeway = width.leeway)
  }

  # Same for indirect
  if ( 'NMRResonance1D' %in% class(indirect.peaks) ) {
    indirect <- indirect.peaks
  } else {
    indirect <- nmrresonance_1d(indirect.peaks, sf = sf[2], id = id,
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


#========================================================================>
#  Lineshape and area calculations
#========================================================================>



#' @rdname f_lineshape
#' @export
setMethod("f_lineshape", "NMRScaffold2D",
  function(object, sf = c(nmroptions$direct$sf, nmroptions$indirect$sf),
           sum.peaks = TRUE, include.id = FALSE, components = 'rr/ii') {

    # Checking to make sure that sweep frequency is defined
    err <- paste('"sf" must be provided as input or set using',
                 'nmroptions$direct$sf and nmroptions$indirect$sf')
    if ( is.null(sf) ) stop(err)

    if ( length(sf) == 1 ) sf <- c(sf, sf)

    # Defining which components to return
    components <- rev(sort(strsplit(a, '[^ri]+', perl = TRUE)[[1]]))

    # Checking 2D components
    err <- paste('"component" argument must consist of two-character codes',
                 'and possibly a separator, e.g., "rr/ii" or "rr ri ir ii"')
    if ( any(! components %in% c('rr', 'ri', 'ir', 'ii')) ) stop(err)

    f_out <- function(y) {
      d <- as_tibble(y)[, components]
      if ( length(components) == 1 ) {
        d[[components]]
      } else if ( length(components == 2) ) { 
        colnames(d) <- c('r', 'i')
        as_cmplx1d(d)
      } else if ( return.i ) {
        as_cmplx2d(d)
      }
    }

    columns <- c('position', 'width', 'height', 'fraction.gauss')
    peaks <- peaks(object)
    parameters <- as.matrix(peaks[, columns])

    direct.peaks <- peaks$dimension == 'direct'
    indirect.peaks <- peaks$dimension == 'indirect'
    i.dim

    # Converting peak width to ppm
    parameters[, 2] <- parameters[, 2]/sf

    # If peaks are to be summed, just feed all parameters into the Rcpp function
    if ( sum.peaks ) {
      out <- function(x) {
        p <- as.vector(t(parameters))
        y <- matrix(0, nrow = length(x), ncol = 2)
        .Call('_rnmrfit_lineshape_1d', PACKAGE = 'rnmrfit', x, y, p)
        f_out(cmplx1(r = y[,1], i = y[,2]))
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
          p <- as.vector(p)
          y <- matrix(0, nrow = length(x), ncol = 2)
          .Call('_rnmrfit_lineshape_1d', PACKAGE = 'rnmrfit', x, y, p)
          f_out(cmplx1(r = y[,1], i = y[,2]))
        }
      })

      # Adding functions as a column
      out$f <- functions
    }

    out
  })
