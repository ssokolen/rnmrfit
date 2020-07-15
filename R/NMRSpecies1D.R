# Definition of a class structure for collection of 1D resonances.



#==============================================================================>
#  NMRSpecies1D -- collection of NMRResonance1D objects
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR species.
#' 
#' Essentially, this class is used to collect and define area relations between
#' multiple resonances as part of a single species.
#' 
#' @slot resonances A list of NMRResonance1D objects.
#' @slot connections A data.frame relating the areas of the resonances.
#' @slot connections.leeway A value specifying how tightly enforced the
#'                          connection constraints on resonance areas should be.
#'                          E.g. a value of = 0 specifies that the area ratios
#'                          are exact, whereas 0.1 specifies that the
#'                          areas of the resonances may differ by +/- 10
#'                          percent from the specified ratios.
#' 
#' @name NMRSpecies1D-class
#' @export
NMRSpecies1D <- setClass("NMRSpecies1D",
  contains = 'NMRScaffold1D',
  slots = c(
    name = 'character',
    id = 'character',
    children = 'list',
    connections = 'data.frame',
    connections.leeway = 'numeric'
  ),
  prototype = prototype(
    name = 'species',
    id = 'species',
    children = list(),
    connections = data.frame(),
    connections.leeway = 0
  )
)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRSpecies1D object
#' 
#' Generates an NMRSpecies1D object from a list of NMRResonance1D objects or a
#' list of character/numeric vectors that can be converted to NMRResonance1D
#' objects. See ?nmrresonance_1d for more details about this conversion.
#' 
#' @param resonances A list of NMRResonance1D objects or a list of
#'                   character/numeric vectors that can be converted to
#'                   NMRResonance1D objects. See ?nmrresonance_1d for more
#'                   details about this conversion. If list elements are named,
#'                   these names will be use to replace resonance ids.
#' @param areas A vector of areas corresponding to the expected area ratios of
#'              the resonances. Set to NULL by default, signifying no fixed
#'              constraints.
#' @param id A string specifying species name. If left empty, a name is
#'           automatically generated from the resonance names.
#' @param connections.leeway A value specifying how tightly enforced the
#'                           connection constraints on resonance areas should
#'                           be. E.g. a value of = 0 specifies that the area
#'                           ratios are exact, whereas 0.1 specifies that the
#'                           areas of the resonances may differ by +/- 10
#'                           percent from the specified ratios.
#' @param ... Options passed to nmrresonance_1d if resonances are being
#'            converted from character/numeric vectors. See ?nmrresonance_1d for
#'            more details.
#' 
#' @return An NMRSpecies1D object.
#' 
#' @export
nmrspecies_1d <- function(resonances, areas = NULL, id = NULL, 
                          connections.leeway = 0, ...) {


  #---------------------------------------
  # Generating list of resonances
  resonances.list <- list()

  # If the original resonances aren't a list, place them into a list
  if ( class(resonances) == 'character' ) resonances <- as.list(resonances)
  else if ( class(resonances) != 'list' ) resonances <- list(resonances)

  for (i in 1:length(resonances)) {

    resonance <- resonances[[i]]

    if ( class(resonance) == 'NMRSpecies1D' ) {
      err <- paste("An NMRSpecies1D can't be constructed from other",
                   "NMRSpecies1D objects. Use resonances() to first extract",
                   "the resonance list before creating a new object.")
      stop(err)
    }
    # If the object is already an NMRResonance1D object, add it directly
    else if ( class(resonance) == 'NMRResonance1D' ) {
      resonances.list <- c(resonances.list, resonance)
    }
    # Otherwise, feed it into the nmrresonance_1d constructor
    else {
      resonances.list <- c(resonances.list, nmrresonance_1d(resonance, ...))
    }
        
    # Modifying id if provided
    resonance.id <- names(resonances)[i]
    if (! is.null(resonance.id) ) resonances.list[[i]]@id <- resonance.id
  }

  #---------------------------------------
  # Defining connections if areas provided

  # Fetching resonance ids from list
  valid.ids <- unlist(lapply(resonances.list, function(o) o@id))

  if (! is.null(areas) ) {
    # Checking that areas correspond to resonances
    err <- paste('Either "areas" vector must have names or the length of',
                 '"areas" vector must match length of resonances.')
    if (! is.null(names(areas)) ) ids <- names(areas)
    else if ( length(areas) == length(valid.ids) ) ids <- valid.ids
    else stop(err) 

    # Checking that area names are valid
    err <- 'Names of "areas" vector must be valid resonance ids.'
    if ( any(! ids %in% valid.ids) ) stop(err)

    # Checking length
    err <- '"areas" vector must be of length 2 or more to add constraints.'
    if ( length(ids) < 2 ) stop(err)

    # Generating connections data frame
    n <- length(areas)
    index.1 <- 1:(n - 1)
    index.2 <- 2:n
    connections <- data.frame(resonance.1 = ids[index.1], 
                              resonance.2 = ids[index.2],
                              area.ratio = areas[index.2]/areas[index.1])
    rownames(connections) <- 1:nrow(connections) 
  } 
  else {
    connections = data.frame()
  }

  # Generating id if it doesn't exist
  if ( is.null(id) ) id <- paste(valid.ids, collapse = '-')

  new('NMRSpecies1D', id = id, children = resonances.list, 
                      connections = connections, 
                      connections.leeway = connections.leeway)
}
