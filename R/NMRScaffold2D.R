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
)



#==============================================================================>
# Access to direct and indirect components
#==============================================================================>



#------------------------------------------------------------------------------
# Direct

#---------------------------------------
#' Get 1D projection of direct dimension 
#' 
#' Generic convenience method to extract all direct dimension components.
#' 
#' @param object An NMRScaffold2D object.
#' 
#' @name direct
setGeneric("direct", 
  function(object) standardGeneric("direct")
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



#------------------------------------------------------------------------------
# Indirect

#---------------------------------------
#' Get 1D projection of indirect dimension 
#' 
#' Generic convenience method to extract all indirect dimension components.
#' 
#' @param object An NMRScaffold2D object.
#' 
#' @name indirect 
setGeneric("indirect", 
  function(object) standardGeneric("indirect")
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
    id <- ifelse( 'id' %in% slotNames(object), sprintf('(%s)', id(object)), '')

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
        cmplx1(r = d$r, i = d$i)
      } else if ( return.i ) {
        cmplx2(rr = d$rr, ri = d$ri, ir = d$ir, ii = d$i )
      }
    }

    columns <- c('position', 'width', 'height', 'fraction.gauss')

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
