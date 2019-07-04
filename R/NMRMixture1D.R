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
#' NMRMixture1D validity test
#'
validNMRMixture1D <- function(object) {

  species <- object@species

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Checking that all species list items are valid
  for ( specie in species ) {
    logic1 <- class(specie) != 'NMRSpecies1D'
    logic2 <- ! validObject(specie)
    if ( logic1 || logic2 ) {
      valid <- FALSE
      new.err <- paste('All elements of "species" list must be valid',
                       'NMRSpecies1D objects.')
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRMixture1D", validNMRMixture1D)



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
#'                elements are named, these names will be use to replace
#'                species ids.
#' @param ... Options passed to nmrspecies_1d if resonances are being
#'            converted from character/numeric vectors. See ?nmrresonance_1d for
#'            more details.
#' 
#' @return An NMRMixture1D object.
#' 
#' @export
nmrmixture_1d <- function(species, ...) {

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
    if (! is.null(specie.id) ) id(species.list[[i]]) <- specie.id
  }

  #---------------------------------------
  # Resulting mixture object
  out <- new('NMRMixture1D', species = species.list)
}



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
# Peaks

#' @rdname peaks
#' @export
setMethod("peaks", "NMRMixture1D", 
  function(object) {
    peaks.list <- lapply(object@species, peaks, include.id = TRUE)
    do.call(rbind, peaks.list)
  })

#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRMixture1D",
  function(object, value) {

    # Check that data.frame has a species column
    err <- 'New peaks data.frame must value a "species" column.'
    if (! 'species' %in% colnames(value) ) stop(err)

    # Check that new species match current species 
    new.names <- unique(value$species)
    old.names <- unlist(lapply(object@species, id))
    logic <- ! new.names %in% old.names
    wrn <- sprintf('The following species are not defined, ignoring: %s',
                   paste(new.names[logic], collapse = ', '))

    if ( any(logic) ) warning(wrn)

    # Splitting up new values and assigning
    new.peaks <- by(value, value$species, function(d) select(d, -species))
    indexes <- which(old.names %in% new.names)

    for ( i in indexes ) {
      specie <- object@species[[i]]
      peaks(specie) <- new.peaks[[old.names[i]]]
      object@species[[i]] <- specie
    }

    validObject(object)
    object 
  })

#' @rdname update_peaks
setMethod("update_peaks", "NMRMixture1D",
  function(object, peaks, exclusion.level = nmroptions$exclusion$level,
           exclusion.notification = nmroptions$exclusion$notification) {

  # Check that columns match before continuing
  current.peaks <- peaks(object)
  err <- '"peaks" columns must match those of current peaks data.frame.'
  if (! all(colnames(peaks) %in% colnames(current.peaks))) stop(err)

  # Check for missing peaks
  current.ids <- apply(current.peaks[, c('species', 'resonance', 'peak')], 1, 
                       paste, collapse = '-')

  new.ids <- apply(peaks[, c('species', 'resonance', 'peak')], 1, 
                   paste, collapse = '-')
  logic <- ! current.ids %in% new.ids

  if ( any(logic) ) {

    msg <- paste('The following peaks were found outside the data range',
                 'and were therefore excluded:\n',
                  paste(current.ids[logic], collapse = ', '))

    # Expanding message based on level
    if ( exclusion.level == 'species' ) {
      removed.species <- unique(current.peaks$species[logic])

      msg <- paste(msg, 
                   '\nBased on the current exclusion.level, the following',
                   'species were further excluded:\n',
                    paste(removed.species, collapse = ', '))

      # Removing resonances from updated peaks
      peaks <- filter(peaks, ! species %in% removed.species)
    }
    else if ( exclusion.level == 'resonance' ) {
      ids <- paste(current.peaks$species, current.peaks$resonance, sep = '-')
      removed.resonances <- unique(ids[logic])

      msg <- paste(msg, 
                   '\nBased on the current exclusion.level, the following',
                   'resonances were further excluded:\n',
                   paste(removed.resonances, collapse = ', '))

      # Removing resonances from updated peaks
      peaks <- filter(peaks, ! resonance %in% removed.resonances)
    }

    # Issue notification as requested
    f.error <- function(x) {
      msg <- paste('"exclusion.notification" must be one "none", "message",',
                   '"warning", or "stop"')
      stop(msg)
    }

    f.notification = switch(exclusion.notification, none = identity,
                            message = message, warning = warning, stop = stop,
                            f.error)
    f.notification(msg)
  }

  # Initializing removal index
  indexes <- c()

  # Updating peaks
  for ( i in 1:length(object@species) ) {
    species <- object@species[[i]]
    id <- species@id
    sub.peaks <- peaks %>% filter(species == id) %>% select(-species)

    if ( nrow(sub.peaks) == 0 ) indexes <- c(indexes, i) 

    species <- update_peaks(species, sub.peaks,
                            exclusion.level = exclusion.level,
                            exclusion.notification = 'none')
    object@species[[i]] <- species
  }

  if ( length(indexes) > 0 ) object@species <- object@species[-indexes]

  object
})


#------------------------------------------------------------------------------
# Couplings

#' @rdname couplings
#' @export
setMethod("couplings", "NMRMixture1D", 
  function(object) {
    couplings.list <- lapply(object@species, couplings, include.id = TRUE)
    do.call(rbind, couplings.list)
  })



#------------------------------------------------------------------------------
# Bounds

#' @rdname bounds
#' @export
setMethod("bounds", "NMRMixture1D", 
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
setMethod("set_general_bounds", "NMRMixture1D",
  function(object, ...) {
    object@species <- lapply(object@species, set_general_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRMixture1D",
  function(object, ...) {
    object@species <- lapply(object@species, set_general_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRMixture1D",
  function(object, ...) {
    object@species <- lapply(object@species, set_general_bounds, ...)
    object
  })

#------------------------------------------------------------------------------
#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRMixture1D",
  function(object, ...) {
    object@species <- lapply(object@species, set_general_bounds, ...)
    object
  })
