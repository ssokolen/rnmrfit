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

  # Updating
  peaks(object) <- peaks
  
  object
})



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

    columns <- c('position', 'width', 'height', 'fraction.gauss')
    peaks <- peaks(object, include.id = include.id)
    parameters <- as.matrix(peaks[, columns])

    columns <- c("mixture", "species", "resonance")
    descriptors <- as.matrix(peaks[, which(colnames(peaks) %in% columns)])
    descriptors <- apply(descriptors, 1, paste, collapse = "-")

    i.res <- as.integer(factor(descriptors)) - 1
    i.res <- rep(i.res, each = 4)

    i.dim <- as.integer(factor(peaks$dimension, 
                               levels = c('direct', 'indirect'))) - 1
    i.dim <- rep(i.dim, each = 4)

    # Converting peak width to ppm
    logic <- peaks$dimension == 'direct'
    parameters[logic, 2] <- parameters[logic, 2]/sf[1]
    parameters[!logic, 2] <- parameters[!logic, 2]/sf[2]
    
    # The overall function is composed of two parts -- the Rust wrapper that
    # calculates values for all dimension and then the R formatter that 
    # selects which of these dimensions to output

    #---------------------------------------
    # First, defining how to format the output 

    components <- rev(sort(strsplit(components, '[^ri]+', perl = TRUE)[[1]]))
    err <- paste('"component" argument must consist of two-character codes',
                 'and possibly a separator, e.g., "rr/ii" or "rr ri ir ii"')
    if ( any(! components %in% c('rr', 'ri', 'ir', 'ii')) ) stop(err)
    
    if ( length(components) == 1 ) {
      f_out <- function(y) {
        d <- as_tibble(y)[, components]
        d[[components]]
      }
    } else if ( length(components) == 2 ) {
      f_out <- function(y) {
        d <- as_tibble(y)[, components]
        cmplx1(r = d[, 1], i = d[, 2])
      }
    } else {
      f_out <- function(y) {
        d <- as_tibble(y)[, components]
        cmplx2(rr = d$rr, ri = d$ri, ir = d$ir, ii = d$ii )
      }
    }

    #---------------------------------------
    # Then, defining wrapper to incorporate the formatting

    f_gen <- function(p, i.res, i.dim) {
      force(p)
      function(x1, x2) {

        n <- as.integer(length(x1))

        y <- .Call("eval_2d_wrapper",        
          x_direct = as.double(x1),
          x_indirect = as.double(x2),
          y = as.double(rep(0, n*4)),
          resonances = as.integer(i.res),
          dimensions = as.integer(i.dim),
          knots = as.double(0),
          p = as.double(as.vector(t(parameters))),
          n = n,
          nl = as.integer(length(parameters)),
          nb = as.integer(0),
          np = as.integer(0),
          nk = as.integer(0)
        )

        f_out(cmplx2(rr = y[1:n], ri = y[(n+1):(2*n)], 
                     ir = y[(2*n+1):(3*n)], ii = y[(3*n+1):(4*n)]))
      }
    }

    #---------------------------------------
    # Finally, output either a single function or a tibble split by resonances

    if ( sum.peaks ) {
      p <- as.vector(t(parameters))
      out <- f_gen(p, i.res, i.dim)
    } 
    # Otherwise, generate a tbl_df data frame
    else {
      columns <- c("mixture", "species", "resonance")
      out <- peaks[, which(colnames(peaks) %in% columns)]
      out <- out[!duplicated(descriptors), ]
      out <- as_tibble(out)

      # The split is by unique resonance rather than row
      par <- split(parameters, i.res)
      res <- split(i.res, i.res)
      dim <- split(i.dim, i.res)

      # Generating a list of functions, each with their parameters enclosed
      functions <- pmap(list(par, res, dim), function (p, r, d) {
        # Split on matrix flattens parameters by column rather than row,
        # so re-formatting it here
        p <- as.vector(t(matrix(p, ncol = 4)))
        f_gen(p, rep(0, length(r)), d)
      })

      # Adding functions as a column
      out$f <- functions
    }

    out
  })



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
