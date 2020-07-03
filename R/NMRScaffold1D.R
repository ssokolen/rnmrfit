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

    columns <- c('position', 'width', 'height', 'fraction.gauss')
    peaks <- peaks(object)
    parameters <- as.matrix(peaks[, columns])

    # Converting peak width to ppm
    parameters[, 2] <- parameters[, 2]/sf

    # The overall function is compose of two parts -- the Rust wrapper that
    # calculates values for all dimension and then the R formatter that 
    # selects which of these dimensions to output

    #---------------------------------------
    # First, defining how to format the output 

    return.r <- grepl('r', tolower(components))
    return.i <- grepl('i', tolower(components))

    err <- '"components" must have at least one of either "r" or "i"'
    if ( return.r && return.i ) f_out <- function(y) {y}
    else if ( return.r ) f_out <- function(y) {Re(y)}
    else if ( return.i ) f_out <- function(y) {Im(y)}
    else stop(err)

    #---------------------------------------
    # Then, defining wrapper to incorporate the formatting

    f_gen <- function(p) {
      force(p)
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

        f_out(cmplx1(r = y[1:n], i = y[(n+1):(2*n)]))
      }
    }

    #---------------------------------------
    # Finally, output either a single function or a tibble split by peaks

    if ( sum.peaks ) {
      p <- as.vector(t(parameters))
      out <- f_gen(p)
    } 
    # Otherwise, generate a tbl_df data frame
    else {
      out <- as_tibble(peaks[, which(! colnames(peaks) %in% columns)])
      parameters <- split(parameters, 1:nrow(parameters))
      
      # Generating a list of functions, each with their parameters enclosed
      functions <- map(parameters, function (p) {
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
    f <- function(g) {
      intensity <- g[[1]](direct.shift) + baseline
      data.frame(direct.shift = direct.shift, 
                 intensity = vec_cast(intensity, complex())) 
    }

    # Note that the unpack/pack functions are used to avoid bind_row errors

    # And apply it for every peak
    d %>% 
      group_by_if(~ ! is.list(.)) %>% 
      do( f(.$f) ) %>%
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
