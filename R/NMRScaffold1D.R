# Definition of a super-class for 1D resonance data.



#==============================================================================>
#  NMRScaffold1D -- super-class for Resonance, Species, Mixture and Fit
#==============================================================================>



#------------------------------------------------------------------------------
#' Super-class for all 1D peak descriptions.
#' 
#' This class is not meant to be used directly. Instead, it provides a set of
#' common methods for NMRResonance1D, NMRSpecies1D, NMRMixture1D, and NMRFit1D
#' objects. As such, it has no slots and all of its methods are meant to be
#' inherited by the objects listed above.
#' 
#' @name NMRScaffold1D-class
#' @export
NMRScaffold1D <- setClass("NMRScaffold1D", 
  contains = "NMRScaffold"
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
  id <- id(object) 
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
  for ( column in columns ) {
    range <- paste('(', lower[, column], ', ', 
                        upper[, column], ')', sep = '')
    peaks[, column] <- range
  }

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
# Validation methods 
#==============================================================================>



#------------------------------------------------------------------------------
#' @rdname conforms
#' @export
setMethod("check_conformity", "NMRScaffold1D",
  function(object, nmrdata, error = TRUE) {

  out <- TRUE

  # Checking nmrdata
  if ( class(nmrdata) != 'NMRData1D' ) {
    err <- '"nmrdata" must be a valid NMRData1D object.'
    if ( error ) stop(err)
    out <- FALSE
  }

  # Checking that data captures all defined peaks
  d <- processed(nmrdata)
  peaks <- peaks(object) 

  logic <- (peaks$position < min(d$direct.shift)) | 
           (peaks$position > max(d$direct.shift))

  if ( any(logic) ) {
    err <- "Some peaks fall outside the data's chemical shift"
    if ( error ) (stop(err))
    out <- FALSE
  }

  out
})



#==============================================================================>
#  Initialization functions (generating parameter estimates based on data)
#==============================================================================>



#------------------------------------------------------------------------------
#' Initialize peak heights of an NMRScaffold1D object
#' 
#' Generates peak height estimates based on spectral data. At this point, there
#' is just one approach: take peak height as the intensity of the data at the
#' current position of the peak. There are plans to develop more sophisticated
#' approaches in the future.
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
  check_conformity(object, nmrdata)

  d <- processed(nmrdata)
  peaks <- peaks(object) 

  # Using loess with very light smoothing
  d <- data.frame(x = d$direct.shift, y = Re(d$intensity))
  m <- loess(y ~ x, data = d, span = 0.1)

  # Generating heights from prediction
  peaks$height <- predict(m, data.frame(x = peaks$position))

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
#' @param direct.sf Sweep frequency (in MHz) -- needed to convert peak widths
#'                  from Hz to ppm. In most cases, it is recommended to set a
#'                  single default value using nmroptions$direct$sf = ..., but
#'                  an override can be provided here.
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



#' @rdname values
#' @export
setMethod("values", "NMRScaffold1D",
  function(object, direct.shift, direct.sf = nmroptions$direct$sf, 
           sum.level = "all", domain = 'r/i', use.cmplx1 = FALSE) {

  # Checking to make sure that sweep frequency is defined
  err <- '"direct.sf" must be provided as input or set using nmroptions$direct$sf'
  if ( is.null(direct.sf) ) stop(err)

  # Generating components to work with a consistent basis
  components <- components(object, sum.level)

  # Function to apply to each component
  data.columns <- c('position', 'width', 'height', 'fraction.gauss')

  f_values <- function(object) {

    # Converting peak width to ppm
    peaks <- peaks(object)
    parameters <- as.matrix(peaks[, data.columns])
    parameters[, 2] <- parameters[, 2]/direct.sf

    p <- as.vector(t(parameters))
    n <- as.integer(length(direct.shift))
    
    y <- .Call(
      "eval_1d_wrapper",        
      x = as.double(direct.shift),
      y = as.double(rep(0, n*2)),
      knots = as.double(0),
      p = as.double(as.vector(p)),
      n = n,
      nl = as.integer(length(p)),
      nb = as.integer(0),
      np = as.integer(0),
      nk = as.integer(0)
    )
    y <- cmplx1(r = y[1:n], i = y[(n+1):(2*n)])

    tibble(direct.shift = direct.shift, intensity = y)
  }

  # And apply it to every component
  d <- components %>% 
    group_by(across(where(~ !is.list(.)))) %>%
    summarise( f_values(component[[1]]) ) %>%
    ungroup()

  # If all components were selected, drop identifiers
  if ( sum.level == "all" ) d <- select(d, direct.shift, intensity)

  # Select output for intensity
  d$intensity <- domain(d$intensity, domain, use.cmplx1)

  d
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
