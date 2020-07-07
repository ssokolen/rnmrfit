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
    children = 'list'
  ),
  prototype = prototype(
    name = 'mixture',
    id = 'mixture',
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
#' @param species A list of NMRSpecies2D objects. If list elements are named,
#'                these names will be use to replace species ids.
#' @param id An optional string specifying mixture name. If left empty, the
#'           mixture name is left as the default "mixture"
#' @param ... Currently ignored.
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
  out <- new('NMRMixture2D', id = id, children = species.list)
}
