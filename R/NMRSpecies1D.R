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
    id = 'character',
    resonances = 'list',
    connections = 'data.frame',
    connections.leeway = 'numeric'
  ),
  prototype = prototype(
    id = 'species',
    resonances = list(),
    connections = data.frame(),
    connections.leeway = 0
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRSpecies1D validity test
#'
validNMRSpecies1D <- function(object) {

  resonances <- object@resonances
  connections <- object@connections 
  connections.leeway <- object@connections.leeway

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Checking name
  if ( length(id) != 1 ) {
    valid <- FALSE
    new.err <- '"name" must be a character vector of length 1.'
    err <- c(err, new.err)
  }

  #---------------------------------------
  # Checking that all resonance list items are valid
  for ( resonance in resonances ) {
    logic1 <- class(resonance) != 'NMRResonance1D'
    logic2 <- ! validObject(resonance)
    if ( logic1 || logic2 ) {
      valid <- FALSE
      new.err <- paste('All elements of "resonances" list must be valid',
                       'NMRResonance1D objects')
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking connections 
  if ( nrow(connections) > 0 ) {

    valid.columns <- c('resonance.1', 'resonance.2', 'area.ratio')
    if (! identical(colnames(connections), valid.columns) ) {
      valid <- FALSE
      new.err <- sprintf('"connections" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      err <- c(err, new.err)
    }

  }

  #---------------------------------------
  # Checking connections leeway
  if ( (connections.leeway < 0) || (connections.leeway >= 1) ) {
    new.err <- '"connections.leeway" must be in the range [0, 1).'
    valid <- FALSE
    err <- c(err, new.err)
  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRSpecies1D", validNMRSpecies1D)



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
#' @param areas A vector of areas corresponding to the expected areas of the
#'              resonances. Set to NULL by default, signifying no fixed
#'              constraings.
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
    if (! is.null(resonance.id) ) id(resonances.list[[i]]) <- resonance.id
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

  new('NMRSpecies1D', id = id, resonances = resonances.list, 
                      connections = connections, 
                      connections.leeway = connections.leeway)
}



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
# Id

#' @rdname id
#' @export
setMethod("id", "NMRSpecies1D", 
  function(object) object@id
  )

#' @rdname id-set
#' @export
setReplaceMethod("id", "NMRSpecies1D",
  function(object, value) {
    object@id <- as.character(value)
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Peaks

#' @rdname peaks
#' @export
setMethod("peaks", "NMRSpecies1D", 
  function(object, include.id = FALSE) {
    peaks.list <- lapply(object@resonances, peaks, include.id = TRUE)
    peaks <- do.call(rbind, peaks.list)
    if ( include.id && (nrow(peaks) > 0) ) cbind(species = object@id, peaks)
    else peaks
  })

#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRSpecies1D",
  function(object, value) {

    # Check that data.frame has a resonance column
    err <- 'New peaks data.frame must value a "resonance" column.'
    if (! 'resonance' %in% colnames(value) ) stop(err)

    # Check that new resonances match current resonances
    new.names <- unique(value$resonance)
    old.names <- unlist(lapply(object@resonances, id))
    logic <- ! new.names %in% old.names
    wrn <- sprintf('The following resonances are not defined, ignoring: %s',
                   paste(new.names[logic], collapse = ', '))

    if ( any(logic) ) warning(wrn)

    # Splitting up new values and assigning
    new.peaks <- by(value, value$resonance, function(d) select(d, -resonance))
    indexes <- which(old.names %in% new.names)

    for ( i in indexes ) {
      resonance <- object@resonances[[i]]
      peaks(resonance) <- new.peaks[[old.names[i]]]
      object@resonances[[i]] <- resonance
    }

    validObject(object)
    object 
  })

#' @rdname update_peaks
setMethod("update_peaks", "NMRSpecies1D",
  function(object, peaks, exclusion.level = nmroptions$exclusion$level,
           exclusion.notification = nmroptions$exclusion$notification) {

  # Check that columns match before continuing
  current.peaks <- peaks(object)
  err <- '"peaks" columns must match those of current peaks data.frame.'
  if (! all(colnames(peaks) %in% colnames(current.peaks))) stop(err)

  # Check for missing peaks
  current.ids <- apply(current.peaks[, c('resonance', 'peak')], 1, 
                       paste, collapse = '-')

  new.ids <- apply(peaks[, c('resonance', 'peak')], 1, 
                   paste, collapse = '-')
  logic <- ! current.ids %in% new.ids

  if ( any(logic) ) {

    msg <- paste('The following peaks were found outside the data range',
                 'and were therefore excluded:\n',
                  paste(current.ids[logic], collapse = ', '))

    # Expanding message based on level
    if ( exclusion.level %in% c('resonance', 'species') ) {
      removed.resonances <- unique(current.peaks$resonance[logic])

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
  for ( i in 1:length(object@resonances) ) {
    resonance <- object@resonances[[i]]
    id <- resonance@id
    sub.peaks <- peaks %>% filter(resonance == id) %>% select(-resonance)

    if ( nrow(sub.peaks) == 0 ) indexes <- c(indexes, i) 

    resonance <- update_peaks(resonance, sub.peaks,
                              exclusion.level = exclusion.level,
                              exclusion.notification = 'none')
    object@resonances[[i]] <- resonance
  }

  if ( length(indexes) > 0 ) object@resonances <- object@resonances[-indexes]

  object
})


#------------------------------------------------------------------------------
# Couplings

#' @rdname couplings
#' @export
setMethod("couplings", "NMRSpecies1D", 
  function(object, include.id = FALSE) {
    couplings.list <- lapply(object@resonances, couplings, include.id = TRUE)
    couplings <- do.call(rbind, couplings.list)
    if ( include.id && (nrow(couplings) > 0) ) {
      cbind(species.1 = object@id, species.2 = object@id, couplings)
    }
    else couplings
  })



#------------------------------------------------------------------------------
# Bounds

#' @rdname bounds
#' @export
setMethod("bounds", "NMRSpecies1D", 
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
setMethod("set_general_bounds", "NMRSpecies1D",
  function(object, ...) {
    object@resonances <- lapply(object@resonances, set_general_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRSpecies1D",
  function(object, ...) {
    object@resonances <- lapply(object@resonances, set_offset_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRSpecies1D",
  function(object, ...) { 
    object@resonances <- lapply(object@resonances, set_conservative_bounds, ...)
    object
  })

#------------------------------------------------------------------------------
#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRSpecies1D",
  function(object, ...) {
    object@resonances <- lapply(object@resonances, set_peak_type, ...)
    object
  })


