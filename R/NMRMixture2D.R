# Definition of a class structure for collection of 2D species.



#==============================================================================>
#  NMRMixture2D -- collection of NMRSpecies2D objects
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of a mixture of species
#' 
#' In the same way that an NMRSpecies2D object is just a collection of
#' NMRResonance2D objects, an NMRMixture2D object is a collection of
#' NMRSpecies2D objects. This class serves primarily for organizational
#' purposes and as the basis of fitting methods.
#' 
#' @slot species A list of NMRSpecies2D objects.
#' 
#' @name NMRMixture2D-class
#' @export
NMRMixture2D <- setClass("NMRMixture2D",
  contains = 'NMRScaffold2D',
  slots = c(
    species = 'list'
  ),
  prototype = prototype(
    species = list()
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRMixture2D validity test
#'
validNMRMixture2D <- function(object) {

  species <- object@species

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Checking that all species list items are valid
  for ( specie in species ) {
    logic1 <- class(specie) != 'NMRSpecies2D'
    logic2 <- ! validObject(specie)
    if ( logic1 || logic2 ) {
      valid <- FALSE
      new.err <- paste('All elements of "species" list must be valid',
                       'NMRSpecies2D objects.')
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRMixture2D", validNMRMixture2D)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRMixture2D object
#' 
#' Generates an NMRMixture2D object from a list of NMRResonance2D/NMRSpecies2D
#' objects.
#' 
#' @param species A list of NMRSpecies2D objects. If list elements are named,
#'                these names will be use to replace species ids.
#' @param ... Currently ignored.
#' 
#' @return An NMRMixture2D object.
#' 
#' @export
nmrmixture_2d <- function(species, ...) {

  #---------------------------------------
  # Generating list of species 

  # If the original species aren't a list, place them into a list
  if ( class(species) == 'character' ) species <- as.list(species)
  else if ( class(species) != 'list' ) species <- list(species)

  species.list <- list()

  for (i in 1:length(species)) {

    specie <- species[[i]]

    # If the object is already an NMRSpecies2D object, add it directly
    if ( class(specie) == 'NMRSpecies2D' ) {
      species.list <- c(species.list, specie)
    }
    # Unlike 1D, other inputs can't be converted
    else {
      err <- '"species" must be a list of NMRSpecies2D objects.'
      stop(err)
    }
    # Modifying id if provided
    specie.id <- names(species)[i]
    if (! is.null(specie.id) ) id(species.list[[i]]) <- specie.id
  }

  #---------------------------------------
  # Resulting mixture object
  out <- new('NMRMixture2D', species = species.list)
}



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>




#------------------------------------------------------------------------------
# Direct and indirect components



#' @rdname direct
#' @export
setMethod("direct", "NMRMixture2D", 
  function(object) {
    direct.list <- lapply(object@species, direct)
    nmrmixture_1d(direct.list)
  })



#' @rdname indirect
#' @export
setMethod("indirect", "NMRMixture2D", 
  function(object) {
    indirect.list <- lapply(object@species, indirect)
    nmrmixture_1d(indirect.list)
  })




#------------------------------------------------------------------------------
# Peaks

#' @rdname peaks
#' @export
setMethod("peaks", "NMRMixture2D", 
  function(object, ...) {
    peaks.list <- lapply(object@species, peaks, include.id = TRUE)
    do.call(rbind, peaks.list)
  })



#------------------------------------------------------------------------------
# Couplings

#' @rdname couplings
#' @export
setMethod("couplings", "NMRMixture2D", 
  function(object) {
    couplings.list <- lapply(object@species, couplings, include.id = TRUE)
    do.call(rbind, couplings.list)
  })



#------------------------------------------------------------------------------
# Bounds

#' @rdname bounds
#' @export
setMethod("bounds", "NMRMixture2D", 
  function(object, include.id = FALSE) {
    f <- function(o, sublist) bounds(o, include.id = TRUE)[[sublist]]
    lower.list <- lapply(object@species, f, sublist = 'lower')
    upper.list <- lapply(object@species, f, sublist = 'upper')

    lower <- do.call(rbind, lower.list)
    upper <- do.call(rbind, upper.list)

    if ( include.id ) {
      if ( nrow(lower) > 0 ) lower <- cbind(species = object@id, lower)
      if ( nrow(upper) > 0 ) upper <- cbind(species = object@id, upper)
    }

    list(lower = lower, upper = upper)
})



#==============================================================================>
#  Bounds
#==============================================================================>



#------------------------------------------------------------------------------
#' @rdname set_general_bounds
#' @export
setMethod("set_general_bounds", "NMRMixture2D",
  function(object, ...) {
    object@species <- lapply(object@species, set_general_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRMixture2D",
  function(object, ...) {
    object@species <- lapply(object@species, set_general_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRMixture2D",
  function(object, ...) {
    object@species <- lapply(object@species, set_general_bounds, ...)
    object
  })

#------------------------------------------------------------------------------
#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRMixture2D",
  function(object, ...) {
    object@species <- lapply(object@species, set_general_bounds, ...)
    object
  })
