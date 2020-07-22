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
    name = 'character',
    id = 'character',
    sf = 'numeric',
    children = 'list'
  ),
  prototype = prototype(
    name = 'mixture',
    id = 'mixture',
    sf = numeric(0),
    children = list()
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRMixture2D object
#' 
#' Generates an NMRMixture2D object from a list of NMRResonance2D/NMRSpecies2D
#' objects.
#' 
#' @param species A list of NMRSpecies2D or NMRResonance2D objects. If list
#'                elements are named, these names will be use to replace species
#'                ids.
#' @param id An optional string specifying mixture name. If left empty, the
#'           mixture name is left as the default "mixture"
#' @param ... Options passed to nmrspecies_2d if NMRResonance2D objects need to
#'            be converted. See ?nmrspecies_2d for more details.
#' 
#' @return An NMRMixture2D object.
#' 
#' @export
nmrmixture_2d <- function(species, id = "mixture", ...) {

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
    else if ( class(specie) == 'NMRResonance2D' ) {
      species.list <- c(species.list, nmrspecies_2d(specie, ...))
    }
    # Unlike 1D, other inputs can't be converted
    else {
      err <- '"species" must consist of NMRSpecies2D or NMRResonance2D objects.'
      stop(err)
    }
    # Modifying id if provided
    specie.id <- names(species)[i]
    if (! is.null(specie.id) ) id(species.list[[i]]) <- specie.id
  }

  # All species must have the same sweep frequency
  sf <- map(species.list, ~ .@sf)
  n <- length(sf)
  if ( n > 1 ) {
    logic <- all(unlist(map2(sf[1:(n-1)], sf[2:n], identical)))
    err <- 'All species used to define a mixture must have the same "sf".'
    if (! logic ) stop(err)
  }
  sf <- sf[[1]]

  #---------------------------------------
  # Resulting mixture object
  out <- new('NMRMixture2D', id = id, sf = sf, children = species.list)
}
