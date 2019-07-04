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
    direct <- data.frame(dimension = "direct", getter(object@direct))
    indirect <- data.frame(dimension = "indirect", getter(object@indirect))
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




