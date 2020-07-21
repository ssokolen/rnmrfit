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



#' @rdname values
#' @export
setMethod("values", "NMRScaffold1D",
  function(object, direct.shift = NULL, sum.level = "all", domain = 'r/i', 
           use.cmplx1 = FALSE) {

  if ( is.null(direct.shift) ) {
    positions <- peaks(object)$position
    direct.shift <- seq(min(positions) - 0.2, max(positions) + 0.2,
                        length.out = 500)
  }

  # Generating components to work with a consistent basis
  components <- components(object, sum.level)

  # Function to apply to each component
  data.columns <- c('position', 'width', 'height', 'fraction.gauss')

  f_values <- function(object) {

    # Converting peak width to ppm
    peaks <- peaks(object)
    parameters <- as.matrix(peaks[, data.columns])
    parameters[, 2] <- parameters[, 2]/object@sf

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
  function(object, sum.peaks = FALSE, include.id = FALSE, components = 'r/i') {

  # Defining area function
  f <- function(position, width, height, fraction.gauss) {
    # If fraction is 0, treat as Lorentz
    if ( fraction.gauss == 0 ) {
      pi*width*height/object@sf
    }
    # If fraction is 1, treat as Gauss
    else if ( fraction.gauss == 1) {
      sqrt(2*pi)*width*height/object@sf
    }
    # Else, proceed as Voigt
    else {
      l.width <- width/object@sf
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



#==============================================================================>
# Plotting  
#==============================================================================>



#------------------------------------------------------------------------------
#' Plot NMRScaffold1D object
#' 
#' Generates an interactive plot object using the plotly package.
#' 
#' If the input is an NMResonance1D, NMRSpecies1D, or NMRMixture1D, the peak
#' lines are simply drawn in red. If the input is an NMRFit1D object, then the
#' output features more components -- the original data is plotted as a black
#' line, the fit is plotted in red, the baseline is plotted in blue, the
#' residual in green. The fit can be plotted as a composite of all the peaks,
#' or individually.
#' 
#' @param x An NMRScaffold1D object.
#' @param domain One of either 'r' or 'i' corresponding to either real or
#'               imaginary data. are displayed in separate subplots.
#' @param direct.shift Used to override default selection of chemical shift
#'                     values.
#' @param sum.level One of either 'all', 'species', 'resonance', 'peak' to
#'                  specify whether all peaks should be summed together.
#' @param add.baseline TRUE to add calculated baseline correction (if it exists)
#'                     to each fit.
#' @param add.phase TRUE to add the calculated phase correction (if it exists)
#'                  to the data.
#' 
#' @return A ggplot2 plot.
#' 
#' @export
plot.NMRScaffold1D <- function(x, domain = 'r', direct.shift = NULL, 
                               sum.level = 'species', add.baseline = TRUE, 
                               add.phase = TRUE) { 

  #---------------------------------------
  # Update a couple of logicals
  if ( add.baseline && ("baseline" %in% slotNames(x)) ) add.baseline <- TRUE
  else add.baseline <- FALSE

  if ( add.phase && ("phase" %in% slotNames(x)) ) add.phase <- TRUE
  else add.phase <- FALSE

  if ( ("nmrdata" %in% slotNames(x)) ) add.data <- TRUE
  else add.data <- FALSE

  #---------------------------------------
  # The fit should come last, so tacking on the extras first 

  p <- NULL

  # If there is raw data, add it, baseline, and residual
  if ( add.data ) {
    d <- x@nmrdata
    direct.shift <- values(d)$direct.shift
    
    if ( add.phase ) d <- add_phase(d, phase(x), degrees = FALSE)

    p <- plot(d, domain = domain, legendgroup = 1, color = "black", 
              name = "Raw data")

    # Total fit
    d.fit <- nmrdata_1d_from_scaffold(x, direct.shift = direct.shift)

    d.total.fit <- add_baseline(d.fit, baseline(x), knots(x))

    d.baseline <- d.total.fit
    d.baseline@processed$intensity <-
      d.total.fit@processed$intensity - d.fit@processed$intensity

    p <- lines(d.baseline, p, domain = domain, legendgroup = 2, color = "blue", 
               name = "Baseline")

    d.residual <- d.total.fit
    d.residual@processed$intensity <- 
      d@processed$intensity - d.total.fit@processed$intensity

    p <- lines(d.residual, p, domain = domain, legendgroup = 3, color = "green", 
               name = "Residual")
  }

  #---------------------------------------
  # Select direct.shift based on data if it exists

  if ( is.null(direct.shift) && add.data ) {
    direct.shift <- values(x@nmrdata)$direct.shift
  }

  #---------------------------------------
  # Then scaffold/fit

  components <- components(x, level = sum.level)$component

  # Initialize plot with the first entry
  d <- components[[1]] %>%
    nmrdata_1d_from_scaffold(direct.shift = direct.shift)

  if ( add.baseline ) d <- add_baseline(d, baseline(x), knots(x))
  
  if ( is.null(p) ) {
    p <- plot(d, domain = domain, legendgroup = 4, color = "red", 
              name = components[[1]]@id)
  } else {
    p <- lines(d, p,  domain = domain, legendgroup = 4, color = "red", 
               name = components[[1]]@id)
  }

  # If there are more components, add them on
  if ( length(components) > 1 ) {
    for ( i in 2:length(components) ) {
      d <- components[[i]] %>%
        nmrdata_1d_from_scaffold(direct.shift = direct.shift)

      if ( add.baseline ) d <- add_baseline(d, baseline(x), knots(x))
        
      p <- lines(d, p, domain = domain, legendgroup = 4, color = "red", 
                 name = components[[i]]@id)
    }
  }


  p
}

setMethod("plot", "NMRScaffold1D", plot.NMRScaffold1D)
