# Definition of a class structure for collection of 2D resonances.



#==============================================================================>
#  NMRSpecies2D -- collection of NMRResonance2D objects
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR species.
#' 
#' Essentially, this class is used to collect and define volume relations
#' between multiple resonances as part of a single species.
#' 
#' @slot children A list of NMRResonance2D objects.
#' @slot connections A data.frame relating the volumes of the resonances.
#' @slot connections.leeway A value specifying how tightly enforced the
#'                          connection constraints on resonance volumes should
#'                          be. E.g. a value of = 0 specifies that the volume
#'                          ratios are exact, whereas 0.1 specifies that the
#'                          volumes of the resonances may differ by +/- 10
#'                          percent from the specified ratios.
#' 
#' @name NMRSpecies2D-class
#' @export
NMRSpecies2D <- setClass("NMRSpecies2D",
  contains = 'NMRScaffold2D',
  slots = c(
    name = 'character',
    id = 'character',
    children = 'list',
    connections = 'data.frame',
    connections.leeway = 'numeric'
  ),
  prototype = prototype(
    name = "species",
    id = "species",
    children = list(),
    connections = data.frame(),
    connections.leeway = 0
  )
)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRSpecies2D object
#' 
#' Generates an NMRSpecies2D object from a list of NMRResonance2D objects.
#' 
#' @param resonances A list of NMRResonance2D objects.
#' @param volumes A vector of volumes corresponding to the expected volume
#'                ratios of the resonances. Set to NULL by default, signifying
#'                no fixed constraints.
#' @param id A string specifying species name. If left empty, a name is
#'           automatically generated from the resonance names.
#' @param connections.leeway A value specifying how tightly enforced the
#'                           connection constraints on resonance volumes should
#'                           be. E.g. a value of = 0 specifies that the volume 
#'                           ratios are exact, whereas 0.1 specifies that the
#'                           volumes of the resonances may differ by +/- 10
#'                           percent from the specified ratios.
#' @param ... Currently ignored. 
#' 
#' @return An NMRSpecies2D object.
#' 
#' @export
nmrspecies_2d <- function(resonances, volumes = NULL, id = NULL, 
                          connections.leeway = 0, ...) {


  #---------------------------------------
  # Generating list of resonances
  resonances.list <- list()

  # If the original resonances aren't a list, place them into a list
  if ( class(resonances) != 'list' ) resonances <- list(resonances)

  for (i in 1:length(resonances)) {

    resonance <- resonances[[i]]

    if ( class(resonance) == 'NMRSpecies2D' ) {
      err <- paste("An NMRSpecies2D can't be constructed from other",
                   "NMRSpecies2D objects. Use resonances() to first extract",
                   "the resonance list before creating a new object.")
      stop(err)
    }
    # If the object is already an NMRResonance2D object, add it directly
    else if ( class(resonance) == 'NMRResonance2D' ) {
      resonances.list <- c(resonances.list, resonance)
    }
    # Unlike 1D, other inputs can't be converted
    else {
      err <- '"resonances" must be a list of NMRResonance2D objects.'
      stop(err)
    }
        
    # Modifying id if provided
    resonance.id <- names(resonances)[i]
    if (! is.null(resonance.id) ) id(resonances.list[[i]]) <- resonance.id
  }

  #---------------------------------------
  # Defining connections if volumes provided

  # Fetching resonance ids from list
  valid.ids <- unlist(lapply(resonances.list, function(o) id(o)))

  if (! is.null(volumes) ) {
    # Checking that volumes correspond to resonances
    err <- paste('Either "volumes" vector must have names or the length of',
                 '"volumes" vector must match length of resonances.')
    if (! is.null(names(volumes)) ) ids <- names(volumes)
    else if ( length(volumes) == length(valid.ids) ) ids <- valid.ids
    else stop(err) 

    # Checking that volume names are valid
    err <- 'Names of "volumes" vector must be valid resonance ids.'
    if ( any(! ids %in% valid.ids) ) stop(err)

    # Checking length
    err <- '"volumes" vector must be of length 2 or more to add constraints.'
    if ( length(ids) < 2 ) stop(err)

    # Generating connections data frame
    n <- length(volumes)
    index.1 <- 1:(n - 1)
    index.2 <- 2:n
    connections <- data.frame(resonance.1 = ids[index.1], 
                              resonance.2 = ids[index.2],
                              volume.ratio = volumes[index.2]/volumes[index.1])
    rownames(connections) <- 1:nrow(connections) 
  } 
  else {
    connections = data.frame()
  }

  # Generating id if it doesn't exist
  if ( is.null(id) ) id <- paste(valid.ids, collapse = '-')

  new('NMRSpecies2D', id = id, children = resonances.list, 
                      connections = connections, 
                      connections.leeway = connections.leeway)
}



#------------------------------------------------------------------------------
# Bounds

#' @rdname bounds
#' @export
setMethod("bounds", "NMRSpecies2D", 
  function(object, include.id = FALSE) {
    f <- function(o, sublist) bounds(o, include.id = TRUE)[[sublist]]
    lower.list <- lapply(object@resonances, f, sublist = 'lower')
    upper.list <- lapply(object@resonances, f, sublist = 'upper')

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
setMethod("set_general_bounds", "NMRSpecies2D",
  function(object, ...) {
    object@resonances <- lapply(object@resonances, set_general_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRSpecies2D",
  function(object, ...) {
    object@resonances <- lapply(object@resonances, set_offset_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRSpecies2D",
  function(object, ...) { 
    object@resonances <- lapply(object@resonances, set_conservative_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRSpecies2D",
  function(object, ...) {
    object@resonances <- lapply(object@resonances, set_peak_type, ...)
    object
  })
