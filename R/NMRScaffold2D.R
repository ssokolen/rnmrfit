# Definition of a super-class for 2D resonance data.



#==============================================================================>
#  NMRScaffold2D -- super-class for Resonance, Species, Mixture and Fit
#==============================================================================>



#------------------------------------------------------------------------------
#' Super-class for all 2D peak descriptions.
#' 
#' This class is not meant to be used directly. Instead, it provides a common
#' framework for methods around visualization and inspection of NMRScaffold2D,
#' NMRSpecies2D, NMRMixture2D, and NMRFit2D objects. Essentially, all of the 2D
#' methods merely wrap around their 1D counterparts, meaning that the same
#' basic approach can be used for all of them.
#'
#' @slot direct An NMRScaffold1D object.
#' @slot indirect An NMRScaffold1D object.
#' 
#' @name NMRScaffold2D-class
#' @export
NMRScaffold2D <- setClass("NMRScaffold2D",
  contains = "VIRTUAL",
  slots = c(
     direct = "NMRScaffold1D",
     indirect = "NMRScaffold1D"
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRScaffold2D validity test
#'
validNMRScaffold2D <- function(object) {

  direct <- object@direct
  indirect <- object@indirect

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Both 1D elements must be valid NMRScaffold1D objects with the same ids
  logic1 <- validObject(direct)
  logic2 <- validObject(indirect)

  if ( logic1 && logic2 ) {

    # If both are valid objects check ids
    if ( direct@id != indirect@id ) {
      valid <- FALSE
      new.err <- '"direct" and "indirect" components must have the same id.'
      err <- c(err, new.err)
    }

  } else {
    
    valid <- FALSE
    new.err <- paste('"direct" and "indirect" components must be valid',
                     'NMRScaffold1D objects.')
    err <- c(err, new.err)
  } 

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRScaffold2D", validNMRScaffold2D)




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

    direct <- object@direct
    indirect <- object@indirect
    id <- direct@id

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
# Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
# Helper functions

#---------------------------------------
#' Combine direct and indirect dimensions
#' 
#' This is an internal function used for all getter functions that output a
#' data.frame object. Essentially, the getter is passed on to the direct and
#' indirect components, a dimension column is added and the resulting objects
#' are stitched together.
#' 
#' @param object An NMRScaffold2D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name combine_getter
setGeneric("combine_getter", 
  function(object, ...) standardGeneric("combine_getter")
)

#' @rdname combine_getter
#' @export
setMethod("combine_getter", "NMRScaffold2D", 
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
#' @param object An NMRScaffold2D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name split_setter
setGeneric("split_setter", 
  function(object, setter, value) standardGeneric("split_setter")
)

#' @rdname split_setter
#' @export
setMethod("split_setter", "NMRScaffold2D", 
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




#------------------------------------------------------------------------------
# Direct

#---------------------------------------
#' Get direct dimension resonance object
#' 
#' Generic convenience method to access the direct dimension resonance
#' component of an NMRScaffold2D object.
#' 
#' @param object An NMRScaffold2D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name direct
#' @export
setGeneric("direct", 
  function(object, ...) standardGeneric("direct")
)

#' @rdname direct
#' @export
setMethod("direct", "NMRScaffold2D", 
  function(object) object@direct
)

#---------------------------------------
#' Set direct dimension resonance object
#' 
#' Generic convenience method to set the direct dimension resonance component
#' of an NMRScaffold2D object.
#' 
#' @param object An NMRScaffold2D object.
#' @param value A valid NMRScaffold1D object.
#' 
#' @name direct-set
#' @export
setGeneric("direct<-", 
  function(object, value) standardGeneric("direct<-")
)

#' @rdname direct-set
#' @export
setReplaceMethod("direct", "NMRScaffold2D",
  function(object, value) {
    object@direct <- value
    validObject(object)
    object 
})



#------------------------------------------------------------------------------
# Indirect

#---------------------------------------
#' Get indirect dimension resonance object
#' 
#' Generic convenience method to access the indirect dimension resonance
#' component of an NMRScaffold2D object.
#' 
#' @param object An NMRScaffold2D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name indirect
#' @export
setGeneric("indirect", 
  function(object, ...) standardGeneric("indirect")
)

#' @rdname indirect
#' @export
setMethod("indirect", "NMRScaffold2D", 
  function(object) object@indirect
)

#---------------------------------------
#' Set indirect dimension resonance object
#' 
#' Generic convenience method to set the indirect dimension resonance component
#' of an NMRScaffold2D object.
#' 
#' @param object An NMRScaffold2D object.
#' @param value A valid NMRScaffold1D object.
#' 
#' @name indirect-set
#' @export
setGeneric("indirect<-", 
  function(object, value) standardGeneric("indirect<-")
)

#' @rdname indirect-set
#' @export
setReplaceMethod("indirect", "NMRScaffold2D",
  function(object, value) {
    object@indirect <- value
    validObject(object)
    object 
})



#------------------------------------------------------------------------------
# Id



#' @rdname id
#' @export
setMethod("id", "NMRScaffold2D", 
  function(object) object@direct@id
)



#' @rdname id-set
#' @export
setReplaceMethod("id", "NMRScaffold2D",
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
setMethod("peaks", "NMRScaffold2D", 
  function(object, include.id = FALSE) {
    combine_getter(object, peaks, include.id = include.id)
})



#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRScaffold2D",
  function(object, value) {
    split_setter(object, `peaks<-`, value)
})



#' @rdname update_peaks
setMethod("update_peaks", "NMRScaffold2D",
  function(object, peaks, exclusion.level = nmroptions$exclusion$level,
           exclusion.notification = nmroptions$exclusion$notification) {
    split_setter(object, update_peaks, peaks, 
                 exclusion.level, exclusion.notification)
})



#------------------------------------------------------------------------------
# Couplings



#' @rdname couplings
#' @export
setMethod("couplings", "NMRScaffold2D", 
  function(object, include.id = FALSE) {
    combine_getter(object, couplings, include.id = include.id)
})



#' @rdname couplings-set
#' @export
setReplaceMethod("couplings", "NMRScaffold2D",
  function(object, value) {
    split_setter(object, `couplings<-`, value)
})



#==============================================================================>
#  Bounds
#==============================================================================>



#' @rdname set_general_bounds
#' @export
setMethod("set_general_bounds", "NMRScaffold2D",
  function(object, position = NULL, height = NULL, width = NULL,
           fraction.gauss = NULL, nmrdata = NULL, widen = FALSE) {

  # Since the nmrdata object is only used to determine maximum and minimum
  # values of chemical shift and intensity, the 2D data is simply split
  # into 1D components
  if (! is.null(nmrdata) ) {

    if ( class(nmrdata) != 'NMRData2D' ) {
      err <- '"nmrdata" must be a valid NMRData2D object.'
      stop(err)
    }
    else {
      validObject(nmrdata)
    }

    columns <- c('direct.shift', 'intensity')
    direct <- new("NMRData1D", processed = processed(nmrdata)[, columns])

    columns <- c('indirect.shift', 'intensity')
    indirect <- new("NMRData1D", processed = processed(nmrdata)[, columns])
  
  } else {
    direct <- NULL
    indirect <- NULL
  }

  object@direct <- set_general_bounds(object@direct, position, height,
                                      width, fraction.gauss, direct, widen)
  object@indirect <- set_general_bounds(object@indirect, position, height,
                                        width, fraction.gauss, indirect, widen)

  validObject(object)
  object
})



#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRScaffold2D",
  function(object, position = NULL, height = NULL, width = NULL, 
           relative = FALSE, widen = FALSE) {

  object@direct <- set_offset_bounds(object@direct, position, height,
                                     width, relative, widen)
  object@indirect <- set_offset_bounds(object@indirect, position, height,
                                       width, relative, widen)

  validobject(object)
  object
})



#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRScaffold2D",
  function(object, position = TRUE,  height = TRUE, width = TRUE, 
           nmrdata = NULL, widen = FALSE) { 

  # The nmrdata object is split up as for general_bounds
  if (! is.null(nmrdata) ) {

    if ( class(nmrdata) != 'NMRData2D' ) {
      err <- '"nmrdata" must be a valid NMRData2D object.'
      stop(err)
    }
    else {
      validObject(nmrdata)
    }

    columns <- c('direct.shift', 'intensity')
    direct <- new("NMRData1D", processed = processed(nmrdata)[, columns])

    columns <- c('indirect.shift', 'intensity')
    indirect <- new("NMRData1D", processed = processed(nmrdata)[, columns])
  
  } else {
    direct <- NULL
    indirect <- NULL
  }

  object@direct <- set_conservative_bounds(object@direct, position, height,
                                           width, direct, widen)
  object@indirect <- set_conservative_bounds(object@indirect, position, height,
                                             width, indirect, widen)

  validObject(object)
  object
})


#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRScaffold2D",
  function(object, peak.type) {

  object@direct <- set_peak_type(object@direct, peak.type)
  object@indirect <- set_peak_type(object@indirect, peak.type)

  validObject(object)
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

    # Defining which components to return
    components <- rev(sort(strsplit(components, '[^ri]+', perl = TRUE)[[1]]))

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
        as_cmplx1(d)
      } else if ( return.i ) {
        as_cmplx2(d)
      }
    }

    columns <- c('position', 'width', 'height', 'fraction.gauss')

    # 2D Resonances must include their ids
    if ( class(object) == 'NMRResonance2D' ) include.id <- TRUE

    peaks <- peaks(object, include.id = include.id)
    parameters <- as.matrix(peaks[, columns])

    # Converting peak width to ppm
    logic <- peaks$dimension == 'direct'
    parameters[logic, 2] <- parameters[logic, 2]/sf[1]
    parameters[!logic, 2] <- parameters[!logic, 2]/sf[2]

    # Defining function generator based on arbitrary subset of parameters
    f_gen <- function(p) {
      force(p)
      function(x1, x2) {
        # Determining unique values
        x.direct <- sort(unique(x1))
        xi.direct <- as.integer(factor(x1), levels = x.direct) - 1

        x.indirect <- sort(unique(x2))
        xi.indirect <- as.integer(factor(x2), levels = x.indirect) - 1

        p <- as.vector(t(parameters))
        y <- matrix(0, nrow = length(x1), ncol = 4)

        i.res <- as.integer(factor(peaks$resonance)) - 1
        i.dim <- as.integer(factor(peaks$dimension, 
                                   levels = c('direct', 'indirect'))) - 1

        .Call('_rnmrfit_lineshape_2d', PACKAGE = 'rnmrfit', 
              x.direct, x.indirect, xi.direct, xi.indirect, y, p, i.res, i.dim)
        f_out(cmplx2(rr = y[,1], ri = y[,2], ir = y[,3], ii = y[,4]))
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



#' @rdname values
#' @export
setMethod("values", "NMRScaffold2D",
  function(object, direct.shift, indirect.shift,
           sf = c(nmroptions$direct$sf, nmroptions$indirect$sf), 
           sum.peaks = TRUE, sum.baseline = FALSE, include.id = FALSE, 
           components = 'rr/ii') {

  err <- '"direct.shift" and "indirect.shift" vectors must be same length'
  if ( length(direct.shift) != length(indirect.shift) ) stop(err)

  # Generating baseline if necessary
  if ( sum.baseline && (class(object) == 'NMRFit2D') ) {
    f <- f_baseline(object, components)
    baseline <- f(direct.shift, indirect.shift)
  } else {
    baseline <- rep(0, length(direct.shift))
  }

  # Output depends on whether peaks are summed or not
  if ( sum.peaks ) {
    # Get function
    f <- f_lineshape(object, sf, sum.peaks, include.id, components)

    # And apply it to specified chemical shifts
    f(direct.shift, indirect.shift) + baseline
  } 
  else {
    # Get data frame of functions
    d <- f_lineshape(object, sf, sum.peaks, include.id, components)

    # Defining function that generates necessary data frame
    f <- function(g) {
      tibble(direct.shift = direct.shift,
             indirect.shift = indirect.shift,
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
#'                  single volume, FALSE to output a data frame of peak area
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

  # Id has to be included at the resonance level
  if ( class(object) == 'NMRResonance2D' ) include.id <- TRUE

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
