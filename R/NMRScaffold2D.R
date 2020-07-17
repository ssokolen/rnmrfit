# Definition of a super-class for 2D resonance data.



#==============================================================================>
#  NMRScaffold2D -- super-class for Resonance, Species, Mixture and Fit
#==============================================================================>



#------------------------------------------------------------------------------
#' Super-class for all 2D peak descriptions.
#' 
#' This class is not meant to be used directly. Instead, it provides a set of
#' common methods for NMRResonance2D, NMRSpecies2D, NMRMixture2D, and NMRFit2D
#' objects. Essentially, all of the 2D methods merely wrap around their 1D
#' counterparts, meaning that the same basic approach can be used for all of
#' them.
#' 
#' @name NMRScaffold2D-class
#' @export
NMRScaffold2D <- setClass("NMRScaffold2D",
  contains = "NMRScaffold",
)



#==============================================================================>
#  Display function
#==============================================================================>



#------------------------------------------------------------------------------
#' Display NMRScaffold2D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRScaffold2D", 
  function(object) {

    direct <- direct(object)
    indirect <- indirect(object)
    id <- id(object)

    cat('\n\n')
    msg <- sprintf('\nAn object of %s class (%s)\n\n', class(object), id)
    cat(strrep('=', nchar(msg) - 3))
    cat(msg)

    cat('------------------------\n')
    cat('In the direct dimension:\n')
    cat('------------------------\n\n')
    show(direct)

    cat('\n\n')
    cat('--------------------------\n')
    cat('In the indirect dimension:\n')
    cat('--------------------------\n\n')
    show(indirect)

})



#==============================================================================>
# Validation methods 
#==============================================================================>



#------------------------------------------------------------------------------
#' @rdname conforms
#' @export
setMethod("check_conformity", "NMRScaffold2D",
  function(object, nmrdata, error = TRUE) {

  out <- TRUE

  # Checking nmrdata
  if ( class(nmrdata) != 'NMRData2D' ) {
    err <- '"nmrdata" must be a valid NMRData2D object.'
    if ( error ) stop(err)
    out <- FALSE
  }

  # Checking that data captures all defined peaks
  d <- processed(nmrdata)
  peaks <- peaks(object) 

  direct.peaks <- filter(peaks, dimension == "direct")
  direct.logic <- (direct.peaks$position < min(d$direct.shift)) | 
                  (direct.peaks$position > max(d$direct.shift)) 


  indirect.peaks <- filter(peaks, dimension == "indirect")
  indirect.logic <- (indirect.peaks$position < min(d$indirect.shift)) | 
                    (indirect.peaks$position > max(d$indirect.shift))

  if ( any(direct.logic) | any(indirect.logic) ) {
    err <- "Some peaks fall outside the data's chemical shift"
    if ( error ) (stop(err))
    out <- FALSE
  }

  out
})



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
#' @param object An NMRScaffold2D object.
#' @param getter Getter function.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name combine_dimensions
setGeneric("combine_dimensions", 
  function(object, ...) standardGeneric("combine_dimensions")
)

#' @rdname combine_dimensions
setMethod("combine_dimensions", "NMRScaffold2D", 
  function(object, getter, ...) {
    out <- map(object@dimensions, getter, ...)
    bind_rows(out, .id = "dimension") 
})



#---------------------------------------
#' Split direct and indirect dimensions
#' 
#' This is an internal function used for all setter functions that input a
#' data.frame object. Essentially, the input value is first split based on a
#' "dimension" column, with the setter being passed on to the direct and
#' indirect components.
#' 
#' @param object An NMRScaffold2D object.
#' @param setter Setter function.
#' @param value Value stored in specified slot.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name split_dimensions
setGeneric("split_dimensions", 
  function(object, setter, value) standardGeneric("split_dimensions")
)

#' @rdname split_dimensions
#' @export
setMethod("split_dimensions", "NMRScaffold2D", 
  function(object, setter, value) {

    # First the input must be a data.frame of some sort
    err <- 'Input value must be a data.frame type object.'
    if (! 'data.frame' %in% class(value) ) stop(err)

    # Second the input must have "dimension" column
    err <- 'Input data.frame must have a "dimension" column.'
    if (! 'dimension' %in% colnames(value) ) stop(err)

    # Third, the dimension column must only contain direct and indirect values
    err <- 'The "dimension" column must contain "direct" and "indirect".'
    entries <- sort(unique(as.character(value$dimension)))
    if (! identical(entries, c('direct', 'indirect')) ) stop(err)

    # If all of the above is met, then split components
    direct <- filter(value, dimension == 'direct') %>% select(-dimension)
    object@dimensions$direct <- setter(object@dimensions$direct, direct)

    indirect <- filter(value, dimension == 'indirect') %>% select(-dimension)
    object@dimensions$indirect <- setter(object@dimensions$indirect, indirect)
    
    validObject(object)
    object
})



#------------------------------------------------------------------------------
# Projection



#' @rdname projection
#' @export
setMethod("projection", "NMRScaffold2D",
  function(object, dimension) {

  slot.names <- slotNames(object)

  # If the object has a dimensions slot, produce contents
  if ( "dimensions" %in% slot.names ) {
    if (! dimension %in% names(object@dimensions) ) {
      err <- sprintf('Dimension "%s" not found', dimension)
      stop(err)
    }

    return(object@dimensions[[dimension]])
  }

  # Otherwise, necessary to construct a new object
  projection <- class(object) %>% str_replace("2", "1")

  # All properties are conserved except children
  f <- function(x) slot(object, x)
  slot.values <- map(slot.names, f) %>%
    set_names(slot.names)

  # nmrdata also has to be projected
  if ( "nmrdata" %in% slot.names ) {
    slot.values$nmrdata <- projection(slot.values$nmrdata, dimension)
  }

  f <- function(...) new(projection, ...)
  projection <- do.call(f, slot.values)

  # The projection is then propogated down to children
  f <- function(x) projection(x, dimension)
  projection@children <- map(projection@children, f) %>%
    set_names(names(projection@children))

  return(projection)
})

#' @rdname projection
#' @export
setMethod("direct", "NMRScaffold2D",
  function(object) {
  
  projection(object, "direct")
})

#' @rdname projection
#' @export
setMethod("indirect", "NMRScaffold2D",
  function(object) {
  
  projection(object, "indirect")
})




#==============================================================================>
#  Initialization functions (generating parameter estimates based on data)
#==============================================================================>



#' @rdname initialize_heights
#' @export
setMethod("initialize_heights", "NMRScaffold2D",
  function(object, nmrdata) {

  # Checking nmrdata
  check_conformity(object, nmrdata)

  # Currently broken
  if ( FALSE ) {
    # Since heights are all dependent anyway, try to get direct dimension
    # height from projection and then leave indirect at default (of 1)
    d <- processed(direct(nmrdata))
    peaks <- peaks(object) 

    # Using loess with very light smoothing
    d <- data.frame(x = d$direct.shift, y = Re(d$intensity))
    m <- loess(y ~ x, data = d, span = 0.1)

    # Generating heights from prediction
    logic <- peaks$dimension == "direct"
    peaks$height[logic] <- predict(m, data.frame(x = peaks$position[logic]))
  }

  d <- processed(nmrdata)
  peaks <- peaks(object) 
  
  logic <- peaks$dimension == "direct"

  y <- Re(d$intensity)
  peaks$height[logic] <- (max(y) + min(y))*0.5

  # Updating
  peaks(object) <- peaks
  
  object
})



#========================================================================>
#  Lineshape and area calculations
#========================================================================>



#' @rdname values
#' @export
setMethod("values", "NMRScaffold2D",
  function(object, direct.shift, indirect.shift,
           direct.sf = nmroptions$direct$sf, 
           indirect.sf = nmroptions$indirect$sf, 
           sum.level = "all", domain = 'rr/ri/ir/ii', use.cmplx1 = FALSE) {

  # Checking to make sure that sweep frequency is defined
  err <- '"direct.sf" must be provided as input or set using nmroptions$direct$sf'
  if ( is.null(direct.sf) ) stop(err)

  err <- '"indirect.sf" must be provided as input or set using nmroptions$indirect$sf'
  if ( is.null(indirect.sf) ) stop(err)

  err <- '"direct.shift" and "indirect.shift" vectors must be same length'
  if ( length(direct.shift) != length(indirect.shift) ) stop(err)

  # Generating components to work with a consistent basis
  components <- components(object, sum.level)

  # Function to apply to each component
  data.columns <- c('position', 'width', 'height', 'fraction.gauss')
  descriptor.columns <- c("mixture", "species", "resonance")

  f_values <- function(object) {

    peaks <- peaks(object)

    # Mapping resonances/dimensions
    logic <- which(colnames(peaks) %in% descriptor.columns)
    descriptors <- as.matrix(peaks[, logic])
    descriptors <- apply(descriptors, 1, paste, collapse = "-")
    dimensions <- peaks$dimension 

    i.res <- as.integer(factor(descriptors))
    i.res <- rep(i.res - 1, each = 4)

    i.dim <- as.integer(factor(dimensions, levels = c('direct', 'indirect')))
    i.dim <- rep(i.dim - 1, each = 4)

    # Converting peak width to ppm
    parameters <- as.matrix(peaks[, data.columns])
    logic <- dimensions == "direct"
    parameters[logic, 2] <- parameters[logic, 2]/direct.sf
    parameters[!logic, 2] <- parameters[!logic, 2]/indirect.sf

    p <- as.vector(t(parameters))
    n <- as.integer(length(direct.shift))
    
    out <- .Call(
      "eval_2d_wrapper",        
      x_direct = as.double(direct.shift),
      x_indirect = as.double(indirect.shift),
      y = as.double(rep(0, n*4)),
      resonances = as.integer(i.res),
      dimensions = as.integer(i.dim),
      knots = as.double(0),
      p = as.double(as.vector(p)),
      n = n,
      nl = as.integer(length(p)),
      nb = as.integer(0),
      np = as.integer(0),
      nk = as.integer(0)
    )

    index.rr <- 1:n
    index.ri <- index.rr + n
    index.ir <- index.ri + n
    index.ii <- index.ir + n

    y <- cmplx2(rr = out[index.rr], ri = out[index.ri],
                ir = out[index.ir], ii = out[index.ii])

    tibble(
      direct.shift = direct.shift, 
      indirect.shift = indirect.shift,
      intensity = y
    )
  }

  # And apply it to every component
  d <- components %>% 
    group_by(across(where(~ !is.list(.)))) %>%
    summarise( f_values(component[[1]]) ) %>%
    ungroup()

  # If all components were selected, drop identifiers
  if ( sum.level == "all" ) {
    d <- select(d, direct.shift, indirect.shift, intensity)
  }

  # Select output for intensity
  d$intensity <- domain(d$intensity, domain, use.cmplx1)

  d
})




#' @rdname areas 
#' @export
setMethod("areas", "NMRScaffold2D", 
          selectMethod("areas", signature = "NMRScaffold1D"))



#------------------------------------------------------------------------
#' Calculate peak volume
#' 
#' Calculate total peak valumes based on peak parameters. Note that the 2D
#' representation consists of multiple "peaks", with at least one in the direct
#' and indirect dimensions. As such, volumes are not calculated per peak, but
#' per resonance.
#' 
#' @param object An NMRScaffold2D object.
#' @param sum.peaks TRUE to add all individual resonances together and output a
#'                  single volume, FALSE to output a data frame of peak volume 
#'                  values.
#' @param include.id TRUE to include id as outer column if outputting data
#'                   frame.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A single overall volume or a data frame of volumes with resaonces
#'         identifier and volume columns.
#' 
#' @name volumes
#' @export
setGeneric("volumes", 
  function(object, ...) standardGeneric("volumes")
)

#' @rdname volumes
#' @export
setMethod("volumes", "NMRScaffold2D",
  function(object, sum.peaks = TRUE, include.id = FALSE) {

  # First, calculate all the individual areas
  areas <- areas(object, sum.peaks = FALSE, include.id = include.id)

  # Volumes are then products of indirect and direct dimension areas
  all.but.height.ids <- colnames(select(areas, -peak, -area))
  all.but.dimension.ids <- colnames(select(areas, -peak, -dimension, -area))

  volumes <- areas %>%
    group_by_at(all.but.height.ids) %>%
    summarize(area = sum(area)) %>%
    group_by_at(all.but.dimension.ids) %>%
    summarize(volume = prod(area))

  # Sum if necessary
  if ( sum.peaks ) sum(volumes$volume)
  else volumes
})



#==============================================================================>
# Plotting  
#==============================================================================>



#------------------------------------------------------------------------------
#' Plot NMRScaffold2D object
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
#' @param x An NMRScaffold2D object.
#' @param domain One of either 'rr', 'ri', 'ir', or 'ii' corresponding to
#'               combinations of real and imaginary data from the direct and
#'               indirect dimensions.
#' @param direct.shift Used to override default selection of chemical shift
#'                     values.
#' @param indirect.shift Similar to direct.shift but for the indirect dimension.
#' @param direct.sf Sweep frequency (in MHz) -- needed to convert peak widths
#'                  from Hz to ppm. In most cases, it is recommended to set a
#'                  single default value using nmroptions$direct$sf = ..., but
#'                  an override can be provided here.
#' @param indirect.sf Similar to direct.sf but for the indirect dimension.
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
plot.NMRScaffold2D <- function(x, domain = 'rr', direct.shift = NULL,
                               indirect.shift = NULL,
                               direct.sf = nmroptions$direct$sf,
                               indirect.sf = nmroptions$indirect$sf,
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
    indirect.shift <- values(d)$indirect.shift
    
    if ( add.phase ) d <- add_phase(d, phase(x), degrees = FALSE)

    p <- plot(d, domain = domain, legendgroup = 1, color = "black", 
              name = "Raw data")

    # Total fit
    d.fit <- nmrdata_2d_from_scaffold(
      x, direct.shift = direct.shift, indirect.shift = indirect.shift,
      direct.sf = direct.sf, indirect.sf = indirect.sf
    )

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
  # Select chemical based on data if it exists

  if ( is.null(direct.shift) && add.data ) {
    direct.shift <- values(x@nmrdata)$direct.shift
  }

  if ( is.null(indirect.shift) && add.data ) {
    direct.shift <- values(x@nmrdata)$indirect.shift
  }

  #---------------------------------------
  # Then scaffold/fit

  components <- components(x, level = sum.level)$component

  # Initialize plot with the first entry
  d <- components[[1]] %>%
    nmrdata_2d_from_scaffold(
      direct.shift = direct.shift, indirect.shift = indirect.shift,
      direct.sf = direct.sf, indirect.sf = indirect.sf
    )

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
        nmrdata_2d_from_scaffold(
          direct.shift = direct.shift, indirect.shift = indirect.shift,
          direct.sf = direct.sf, indirect.sf = indirect.sf
        )

      if ( add.baseline ) d <- add_baseline(d, baseline(x), knots(x))
        
      p <- lines(d, p, domain = domain, legendgroup = 4, color = "red", 
                 name = components[[i]]@id)
    }
  }

  p
}

setMethod("plot", "NMRScaffold2D", plot.NMRScaffold2D)
