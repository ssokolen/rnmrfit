# Definition of a class structure for collection of 1D species.



#==============================================================================>
#  NMRMixture1D -- collection of NMRSpecies1D objects
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of a mixture of species
#' 
#' In the same way that an NMRSpecies1D object is just a collection of
#' NMRResonance1D objects, an NMRMixture1D object is a collection of
#' NMRSpecies1D objects. This class serves primarily for organizational
#' purposes and as the basis of fitting methods.
#' 
#' @slot species A list of NMRSpecies1D objects.
#' 
#' @name NMRMixture1D-class
#' @export
NMRMixture1D <- setClass("NMRMixture1D",
  contains = 'NMRScaffold1D',
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
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRMixture1D object
#' 
#' Generates an NMRMixture1D object from a list of NMRResonance1D/NMRSpecies1D
#' objects or a list of character/numeric vectors that can be converted to
#' NMRResonance1D objects. See ?nmrresonance_1d for more details about this
#' conversion.
#' 
#' @param species A list of NMRSpecies1D objects or other objects that can be
#'                converted to NMRSpecies1D objects. See ?nmrresonance_1d and
#'                ?nmrspecies_1d for more details about this conversion. If list
#'                elements are named, these names will be use to replace species
#'                ids.
#' @param id An optional string specifying mixture name. If left empty, the
#'           mixture name is left as the default "mixture"
#' @param ... Options passed to nmrspecies_1d if resonances are being converted
#'            from character/numeric vectors. See ?nmrresonance_1d for more
#'            details.
#' 
#' @return An NMRMixture1D object.
#' 
#' @export
nmrmixture_1d <- function(species, id = "mixture", ...) {

  #---------------------------------------
  # Generating list of species 

  # If the original species aren't a list, place them into a list
  if ( class(species) == 'character' ) species <- as.list(species)
  else if ( class(species) != 'list' ) species <- list(species)

  species.list <- list()

  for (i in 1:length(species)) {

    specie <- species[[i]]

    # If the object is already an NMRSpecies1D object, add it directly
    if ( class(specie) == 'NMRSpecies1D' ) {
      species.list <- c(species.list, specie)
    }
    # Otherwise, feed it into the nmrspecies_1d constructor
    else {
      species.list <- c(species.list, nmrspecies_1d(specie, ...))
    }
        
    # Modifying id if provided
    specie.id <- names(species)[i]
    if (! is.null(specie.id) ) species.list[[i]]@id <- specie.id
  }

  #---------------------------------------
  # Resulting mixture object
  out <- new('NMRMixture1D', id = id, children = species.list)
}
