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
#' @slot resonances A list of NMRResonance2D objects.
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
#' NMRSpecies2D validity test
#'
validNMRSpecies2D <- function(object) {

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
    logic1 <- class(resonance) != 'NMRResonance2D'
    logic2 <- ! validObject(resonance)
    if ( logic1 || logic2 ) {
      valid <- FALSE
      new.err <- paste('All elements of "resonances" list must be valid',
                       'NMRResonance2D objects')
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking connections 
  if ( nrow(connections) > 0 ) {

    valid.columns <- c('resonance.1', 'resonance.2', 'volume.ratio')
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
setValidity("NMRSpecies2D", validNMRSpecies2D)



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

  new('NMRSpecies2D', id = id, resonances = resonances.list, 
                      connections = connections, 
                      connections.leeway = connections.leeway)
}



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
# Direct and indirect components



#' @rdname direct
#' @export
setMethod("direct", "NMRSpecies2D", 
  function(object) {
    direct.list <- lapply(object@resonances, direct)
    out <- nmrspecies_1d(direct.list, id = id(object))
    
    if ( nrow(object@connections) > 0 ) {
      connections <- object@connections %>%
        mutate(area.ratio = sqrt(volume.ratio)) %>%
        select(-volume.ratio)
      out@connections <- connections
    }

    out
  })



#' @rdname indirect
#' @export
setMethod("indirect", "NMRSpecies2D", 
  function(object) {
    indirect.list <- lapply(object@resonances, indirect)
    out <- nmrspecies_1d(indirect.list, id = id(object))

    if ( nrow(object@connections) > 0 ) {
      connections <-object@connections %>%
        mutate(area.ratio = sqrt(volume.ratio)) %>%
        select(-volume.ratio)
      out@connections <- connections
    }

    out
  })



#------------------------------------------------------------------------------
# Id



#' @rdname id
#' @export
setMethod("id", "NMRSpecies2D", 
  function(object) object@id
  ) 



#' @rdname id-set
#' @export
setReplaceMethod("id", "NMRSpecies2D",
  function(object, value) {
    object@id <- as.character(value)
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Peaks



#' @rdname peaks
#' @export
setMethod("peaks", "NMRSpecies2D", 
  function(object, include.id = FALSE) {
    peaks.list <- lapply(object@resonances, peaks, include.id = TRUE)
    peaks <- do.call(rbind, peaks.list)
    if ( include.id && (nrow(peaks) > 0) ) cbind(species = object@id, peaks)
    else peaks
  })



#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRSpecies2D",
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
setMethod("update_peaks", "NMRSpecies2D",
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
    id <- id(resonance) 
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
setMethod("couplings", "NMRSpecies2D", 
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
