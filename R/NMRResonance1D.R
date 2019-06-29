# Definition of a class structure for 1D resonance data.



#==============================================================================>
#  NMRResonance1D -- peak description 
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR resonance.
#' 
#' Essentially, this class is used to define coupling relationships to group
#' individual peaks into resonances. A resonances is a series of peaks that are
#' constrained by a set of relations between their position, width, and,
#' fraction.gauss parameters as well as overall areas.
#' 
#' @slot id A name to be used for the object in output data.
#' @slot peaks A data.frame describing a series of singlets, with one row per
#'             peak. All peaks are characterized by a position (in ppm), height
#'             (in relative intensity units), width (in ppm or Hz), and
#'             fraction.guass (in percent).
#' @slot couplings A data.frame relating the position and area parameters of the
#'                 peaks, effectively combining singlets into multiplets.
#' @slot couplings.leeway A list specifying how tightly enforced the coupling
#'                        constraints on peak positions, areas, widths, and
#'                        fraction.gauss parameters should be. E.g. position = 0
#'                        specifies that the j coupling constant is exact,
#'                        whereas width = 0.1 specifies that the widths of
#'                        individual peaks may differ by +/- 10 percent.
#' @slot bounds A list of lower and upper bounds on the peak parameters where
#'              both bounds take the same shape as the peaks data.frame.
#' 
#' @name NMRResonance1D-class
#' @export
NMRResonance1D <- setClass("NMRResonance1D",
  contains = 'NMRScaffold1D',
  slots = c(
    id = 'character',
    peaks = 'data.frame',
    couplings = 'data.frame',
    couplings.leeway = 'list',
    bounds = 'list'
  ),
  prototype = prototype(
    id = 'resonance',
    couplings = data.frame(),
    couplings.leeway = list(position = 0, width = 0, 
                            fraction.gauss = 0, area = 0),
    bounds = list(lower = NULL, upper = NULL)
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRResonance1D validity test
#'
validNMRResonance1D <- function(object) {

  id <- object@id
  peaks <- object@peaks
  couplings <- object@couplings
  couplings.leeway <- object@couplings.leeway
  bounds <- object@bounds

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
  # Checking peak column names
  valid.columns <- c('peak', 'position', 'width', 'height', 'fraction.gauss')

  if (! identical(colnames(peaks), valid.columns) ) {
    valid <- FALSE
    new.err <- sprintf('"peaks" must have the following columns: %s',
                       paste(valid.columns, collapse = ', '))
    err <- c(err, new.err)
  }

  #---------------------------------------
  # Checking that lower bounds match peaks
  if (! is.null(bounds$lower) ) {

    logic <- identical(colnames(bounds$lower), valid.columns)
    if (! logic ) {
      valid <- FALSE
      new.err <- sprintf('"bounds$lower" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      err <- c(err, new.err)
    }

    logic1 <- identical(bounds$lower$resonance, peaks$resonance)
    logic2 <- identical(bounds$lower$peak, peaks$peak)
    if (! (logic1 && logic2) ) {
      valid <- FALSE
      new.err <- '"bounds$lower" resonance and peak columns must match "peaks"'
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking that upper bounds match peaks
  if (! is.null(bounds$upper) ) {

    logic <- identical(colnames(bounds$upper), valid.columns)
    if (! logic ) {
      valid <- FALSE
      new.err <- sprintf('"bounds$upper" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      err <- c(err, new.err)
    }

    logic1 <- identical(bounds$upper$resonance, peaks$resonance)
    logic2 <- identical(bounds$upper$peak, peaks$peak)
    if (! (logic1 && logic2) ) {
      valid <- FALSE
      new.err <- '"bounds$upper" resonance and peak columns must match "peaks"'
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking couplings
  if ( nrow(couplings) > 0 ) {

    valid.columns <- c('peak.1', 'peak.2', 'position.difference', 'area.ratio')
    if (! identical(colnames(couplings), valid.columns) ) {
      valid <- FALSE
      new.err <- sprintf('"couplings" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking area and width leeways
  for ( leeway in names(couplings.leeway) ) {
    if ( (couplings.leeway[leeway] < 0) || (couplings.leeway[leeway] >= 1) ) {
      new.err <- sprintf('"couplings.leeway$%s" must be in the range [0, 1).',
                         leeway)
      valid <- FALSE
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRResonance1D", validNMRResonance1D)



#==============================================================================>
#  Helper functions for constructor
#==============================================================================>



#------------------------------------------------------------------------------
#' Parse peak coupling strings
#' 
#' Converts a coupling specification of the form '3.0 d 1.0' into the chemical
#' shift, number of peaks involved and the coupling between them. Currently,
#' the supported codes include s, d, t, q, pentet, sextet, septet, octet, and
#' nonet or any combination of the above. Essentially, the function splits the
#' coupling into a position (number), coupling description (text), and coupling
#' constant (number), parsing the latter two to extract the required codes and
#' numbers. The actual parsing is quite flexible so "1.0 dtpnt 1.5/2/1" will be
#' parsed the same as "1.0 d t pentet 1.5 2 1".
#' 
#' @param coupling.string A character vector with elements of the form '3.0 d
#'                        1.0' or '2 dt 1.0 1.2'
#' 
#' @return A list of three vectors -- chemical.shift, peak numbers, and coupling
#'         constants respectively.
#' 
#' @export
parse_peaks_1d <- function(coupling.string) {

  # Saving original string for later
  original.string <- coupling.string

  # Ensure lowercase
  coupling.string <- tolower(coupling.string)

  # Remove any and all brackets
  coupling.string <- str_replace_all(coupling.string, '[(){}]|\\[\\]', ' ')

  # Replace all possible separators with a single whitespace
  coupling.string <- str_replace_all(coupling.string, '[ ,;#$%_/]+', ' ')

  # Remove whitespace at beginning or end
  coupling.string <- str_trim(coupling.string, 'both')

  # First split directly after first number
  split <- str_split(coupling.string, '(?<=^[0-9.]{1,20})(?![0-9.])', 2)[[1]]
  direct.shift <- split[1]
  direct.shift <- as.numeric(direct.shift)
  coupling.string <- str_trim(split[2], 'both')

  # Then split directly before first remaining number
  split <- str_split(coupling.string, '(?<![0-9.])(?=[0-9])', 2)[[1]]
  codes <- str_trim(split[1], 'both')
  constants <- split[2]

  #---------------------------------------
  # Parsing codes

  # Saving original code text for later
  original <- codes

  # Try every iteration of long names down to three characters
  patterns <- c('pent', 'pnt', 'pentet', 'qui', 'qnt', 'quint', 'quintet',
                'sxt', 'sext', 'sextet', 'spt', 'sept', 'septet', 'hpt',
                'hept', 'heptet', 'oct', 'octet', 'non', 'nonet')
  numbers <- rep(5:9, c(7, 3, 6, 2, 2))

  # Arrange in order of longest to shortest
  index <- order(nchar(patterns), decreasing = TRUE)
  patterns <- patterns[index]
  numbers <- numbers[index]

  replacement <- as.character(numbers)
  names(replacement) <- patterns

  codes <- str_replace_all(codes, replacement)

  # Then parsing single character names
  replacement <- as.character(1:4)
  names(replacement) <- c('s', 'd', 't', 'q')

  codes <- str_replace_all(codes, replacement)

  # Finally, remove any spaces and split by single numbers
  codes <- str_replace_all(codes, '\\s', '')
  codes <- suppressWarnings(as.numeric(str_split(codes, '')[[1]]))

  err <- 'The coupling definition "%s" could not be parsed.'
  if ( any( is.na(codes) ) ) stop(sprintf(err, original))

  # If there is only one singlet, no need to parse constants
  if ( identical(codes, 1) ) {
    return(list(direct.shift = direct.shift, numbers = codes, constants = NA))
  }

  #---------------------------------------
  # Parsing constants

  # Saving original constants text for later
  original <- constants

  # Split by single spaces (previous conversions should have removed issues)
  constants <- suppressWarnings(as.numeric(str_split(constants, ' ')[[1]]))

  err <- 'The coupling constant definition "%s" could not be parsed.'
  if ( any( is.na(constants) ) ) stop(sprintf(err, original))

  # Making sure that the number of codes matches constants
  err <- '"%s" was not parsed with an equal number of couplings and constants.'
  logic <- sum(codes != 1) != sum(! is.na(constants))
  if ( logic ) stop(sprintf(err, original.string))

  # Filling in NA values in case there are singlets
  for (i in which(codes == 1)) {
    if (! is.na(constants[i]) ) {
      constants <- append(constants, NA, i-1)
    }
  }

  list(direct.shift = direct.shift, numbers = codes, constants = constants)
}



#------------------------------------------------------------------------------
#' Split peaks data.frame according to specified coupling.
#' 
#' @param peaks A peaks data.frame from NMRResonance1D object.
#' @param number The number of output peaks per input peak.
#' @param constant The coupling constant of the split.
#' 
#' @return A modified data.frame suitable for NMRResonance1D object.
#' 
#' @export
split_peaks_1d <- function(peaks, number, constant) {

  # Singlets do not require splitting
  if ( number == 1 ) return(peaks)

  # Calculating offsets
  offsets <- constant*(number - 1) * seq(-0.5, 0.5, length = number)

  # Calculating height ratios based on pascal's triangle
  ratios <- choose(number - 1, 0:(number - 1))
  heights <- ratios/sum(ratios)

  # Replicating heights and offsets to match data frame
  n <- nrow(peaks)
  offsets <- rep(offsets, n)
  heights <- rep(heights, n)

  # Replicating each row to match number of output peaks
  peaks <- peaks[rep(1:n, each = number), ]
  peaks$height <- peaks$height*heights
  peaks$position <- peaks$position + offsets

  # Resetting peak number and row names
  peaks <- peaks[order(peaks$position), ]
  peaks$peak <- 1:nrow(peaks)
  rownames(peaks) <- peaks$peak

  peaks
}



#------------------------------------------------------------------------------
#' Enforce coupling relations between peaks.
#' 
#' Adds a set of constraints between specified peaks that fixes the differences
#' between peak positions and the ratios between peak areas to whatever the
#' current positions and peak heights are. These constraints can be relaxed by
#' setting leeway parameters that convert hard equality constraints to soft
#' inequality constraints around a fraction of the fixed values.
#' 
#' @param nmrresonance An NMRResonance1D object to be modified.
#' @param peaks Vector of peak numbers to enforce coupling. The default is to
#'              select all peaks.
#' 
#' @return An NMRResonance1D object with a new set of coupling constraints.
#' 
#' @export
enforce_couplings_1d <- function(nmrresonance, peaks = NULL) {

  # Checking validity
  if ( class(nmrresonance) != 'NMRResonance1D' ) {
    err <- '"nmrresonance" must be a valid NMRResonance1D object.'
    stop(err)
  }
  else {
    validObject(nmrresonance)
  }

  # Filtering peaks
  d <- nmrresonance@peaks
  if ( is.null(peaks) ) peaks <- d$peak
  d <- d[d$peak %in% peaks, ]

  # There must be more than one peak to add couplings
  if ( nrow(d) <= 1 ) return(nmrresonance)

  # Calculating
  index.1 <- 1:(nrow(d) - 1)
  index.2 <- 2:nrow(d)

  differences <- d$position[index.2] - d$position[index.1]
  ratios <- d$height[index.2]/d$height[index.1]
  couplings <- data.frame(peak.1 = d$peak[index.1], peak.2 = d$peak[index.2],
                          position.difference = differences,
                          area.ratio = ratios)

  nmrresonance@couplings <- rbind(nmrresonance@couplings, couplings)
  nmrresonance
}



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRResonance1D object based on simplified peak list
#' 
#' Generates an NMRResonance1D object by converting multiplet definitions into
#' a set of singlets related by constraints on their position and area. Note
#' that peak height parameters are arbitrary at this point. Initial guesses for
#' peak height can be generated during the fit process or overriden manually.
#' 
#' @param peaks A numeric vector of singlet chemical shifts or a character
#'              string specifying multiplets of the form "3 d 1.2". See
#'              ?parse_peaks_1d for more information.
#' @param sf Sweep frequency (in MHz) -- needed to convert coupling constants
#'           from Hz to ppm. In most cases, it is recommended to set a single
#'           default value using nmrsession_1d(sf = ...), but an override can be
#'           provided here.
#' @param id A string specifying resonance name. If left empty, a name is
#'           automatically generated from the peaks argument.
#' @param width Initial estimate of peak width (in Hz). For Voigt lineshapes,
#'              this value is taken as the Lorentzian component, with the
#'              Gaussian component calculated from
#'              peak.width*frac.guass/(1-frac.gauss).
#' @param fraction.gauss Fraction of overall peak width that corresponds to a
#'                       Gaussian lineshape. A value of 0 corresponds to a
#'                       Lorentz peak whereas a value of 1 corresponds to a
#'                       Gaussian peak. Values in between 0 and 1 are modelled
#'                       as a Voigt lineshape but the specific value of
#'                       frac.gauss does not have a physical interpretation.
#' @param position.leeway A fraction specifying how tightly enforced the
#'                        coupling constraints on peak positions, should be.
#'                        E.g. coupling.leeway = 0 specifies that the j coupling
#'                        constant is exact, whereas couping.leeway = 0.1
#'                        specifies that the coupling constant may differ by +/-
#'                        10 percent.
#' @param width.leeway Similar to position.leeway but for peak widths.
#'                     Determines how strictly equal peak widths for all coupled
#'                     peaks are enforced.
#' @param fraction.gauss.leeway Similar to position.leeway but for
#'                              fraction.gauss parameters. Determines how
#'                              strictly equal gauss components are enforced for
#'                              all coupled peaks. A leeway of 0.1, for example,
#'                              would enable one peak to be 24% gauss while
#'                              another peak within the same resonance is 34%.
#' @param area.leeway Similar to position.leeway but for peak areas. Determines
#'                    how strictly the coupling area ratios are enforced.
#' 
#' @return An NMRResonance1D object.
#' 
#' @export
nmrresonance_1d <- function(peaks, sf = nmrsession_1d('sf'), id = NULL, 
                            width = 1, fraction.gauss = 0, 
                            position.leeway = 0, width.leeway = 0, 
                            fraction.gauss.leeway = 0, area.leeway = 0) {

  #---------------------------------------
  # Building peak list

  # Couplings are not added in every case
  add.couplings <- FALSE

  # If peaks is a character, parse coupling information
  if ( is.character(peaks) ) {
    coupling <- parse_peaks_1d(peaks)
    add.constraints <- TRUE

    # Initializing singlet at chemical shift
    if ( is.null(id) ) id <- peaks
    peaks <- data.frame(peak = 1, position = coupling$direct.shift,
                        width = width, height = 1, 
                        fraction.gauss = fraction.gauss)

    # If there is splitting to do, convert constants from Hz to ppm
    if ( any(coupling$number > 1) ) {

      # Checking to make sure that sweep frequency is defined
      err <- '"sf" must be provided as input or set using nmrsession_1d()'
      if ( is.null(sf) ) stop(err)

      # Converting coupling constant from Hz to ppm
      coupling$constant <- coupling$constant/sf

    }

    # Looping through the coupling to split the specified peaks
    for ( i in 1:length(coupling$number) ) {
      peaks <- split_peaks_1d(peaks, coupling$number[i], coupling$constant[i])
    }

    # Set flag to add couplings later
    add.couplings <- TRUE
  }
  # Otherwise, build peaks directly from singlets
  else {
    add.constraints <- FALSE
    middle <- (max(peaks) + min(peaks))/2 
    range <- paste(min(peaks), '..', max(peaks), sep = '')
    if ( is.null(id) ) id <- paste(middle, 'm', range)
    peaks <- data.frame(peak = 1:length(peaks), 
                        position = peaks, width = width, height = 1, 
                        fraction.gauss = fraction.gauss)
  }

  #---------------------------------------
  # Adding coupling definitions

  # Starting with blanks first
  couplings <- data.frame()
  couplings.leeway = list(position = position.leeway, width = width.leeway,
                          fraction.gauss.leeway = fraction.gauss.leeway, 
                          area = area.leeway)

  nmrresonance = new('NMRResonance1D', id = id, peaks = peaks, 
                                       couplings = couplings, 
                                       couplings.leeway = couplings.leeway)

  # And then updating if necessary
  if ( add.couplings ) enforce_couplings_1d(nmrresonance)
  else nmrresonance

}



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>



#' @rdname id
#' @export
setMethod("id", "NMRResonance1D", 
  function(object) object@id
)



#' @rdname id-set
#' @export
setReplaceMethod("id", "NMRResonance1D",
  function(object, value) {
    object@id <- as.character(value)
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Peaks



#' @rdname peaks
#' @export
setMethod("peaks", "NMRResonance1D", 
  function(object, include.id = FALSE) {
    peaks <- object@peaks
    if ( include.id && (nrow(peaks) > 0) ) {
      cbind(resonance = object@id, peaks)
    }
    else peaks
  })



#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRResonance1D",
  function(object, value) {
    object@peaks <- value
    validObject(object)
    object 
  })



#' @rdname update_peaks
setMethod("update_peaks", "NMRResonance1D",
  function(object, peaks, exclusion.level = nmrsession_1d$exclusion$level,
           exclusion.notification = nmrsession_1d$exclusion$notification) {

  # Check that columns match before continuing
  current.peaks <- peaks(object)
  err <- '"peaks" columns must match those of current peaks data.frame.'
  if (! all(colnames(peaks) %in% colnames(current.peaks))) stop(err)

  # Check for missing peaks
  if ( exclusion.level %in% c('resonance', 'species') ) {
    logic <- rep(TRUE, nrow(current.peaks) )
  }
  else {
    current.ids <- current.peaks$peak
    new.ids <- peaks$peak
    logic <- ! current.ids %in% new.ids
  }

  if ( any(logic) ) {

    msg <- paste('The following peaks were found outside the data range',
                 'and were therefore excluded:\n',
                  paste(current.ids[logic], collapse = ', '))

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

  # Setting peaks
  object@peaks <- peaks

  # Updating bounds to make sure they still relate to peaks
  bounds <- object@bounds
  bounds$lower <- bounds$lower[!logic, ]
  bounds$upper <- bounds$upper[!logic, ]
  object@bounds <- bounds

  object
})



#------------------------------------------------------------------------------
# Couplings



#' @rdname couplings
#' @export
setMethod("couplings", "NMRResonance1D", 
  function(object, include.id = FALSE) {
    couplings <- object@couplings
    if ( include.id && (nrow(couplings) > 0) ) {
      cbind(resonance.1 = object@id, resonance.2 = object@id, couplings)
    }
    else couplings
  })



#' @rdname couplings-set
#' @export
setReplaceMethod("couplings", "NMRResonance1D",
  function(object, value) {
    object@couplings <- value
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Bounds

#---------------------------------------
#' Initialize empty set of bounds
#' 
#' Initalize lower and upper bounds with the correct dimensions, but set to
#' -Inf and +Inf respectively
#' 
#' @param object NMRResonance1D object
#' @param overwrite TRUE to overwrite existing bounds (to reset them), FALSE to
#'                  quietly ignore any existing bounds.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return Modified NMRResonance1D object.
#' @name initialize_bounds
setGeneric(".initialize_bounds", 
 function(object,  ...) {
   standardGeneric(".initialize_bounds")
 })

#' @rdname initialize_bounds
setMethod(".initialize_bounds", "NMRResonance1D", 
  function(object, overwrite = FALSE) {

    # Handling bounds if they exist
    bounds <- list(lower = object@bounds$lower, 
                   upper = object@bounds$upper)

    # Selecting default values
    values <- list(lower = -Inf, upper = +Inf)
    columns <- c('position', 'width', 'height', 'fraction.gauss')

    for ( name in names(bounds) ) {
      if ( is.null(bounds[[name]]) || overwrite ) {

        peaks <- object@peaks
        peaks[ , columns] <- values[[name]]

        bounds[[name]] <- peaks
      }
    }

    object@bounds$lower <- bounds$lower
    object@bounds$upper <- bounds$upper

    object
  })



#' @rdname bounds
#' @export
setMethod("bounds", "NMRResonance1D", 
  function(object, include.id = FALSE) {
    bounds <- .initialize_bounds(object)@bounds

    lower <- bounds$lower
    upper <- bounds$upper

    if ( include.id ) {
      if ( nrow(lower) > 0 ) {
        bounds$lower <- cbind(resonance = object@id, lower)
      }

      if ( nrow(upper) > 0 ) {
        bounds$upper <- cbind(resonance = object@id, upper)
      }
    }
    bounds
  })



#' @rdname bounds-set
#' @export
setReplaceMethod("bounds", "NMRResonance1D",
  function(object, value) {
    object@bounds <- value
    validObject(object)
    object 
  })



#==============================================================================>
#  Bounds
#==============================================================================>



#' @rdname set_general_bounds
#' @export
setMethod("set_general_bounds", "NMRResonance1D",
  function(object, position = NULL, height = NULL, width = NULL,
           fraction.gauss = NULL, nmrdata = NULL, widen = FALSE) {
  
  # Initializing bounds
  object <- .initialize_bounds(object)
  lower <- object@bounds$lower
  upper <- object@bounds$upper

  #---------------------------------------
  # Scaling all bounds if nmrdata has been provided
  if (! is.null(nmrdata) ) {

    if ( class(nmrdata) != 'NMRData1D' ) {
      err <- '"nmrdata" must be a valid NMRData1D object.'
      stop(err)
    }
    else {
      validObject(nmrdata)
    }

    processed <- nmrdata@processed
    y.range <- max(Re(processed$intensity)) - min(Re(processed$intensity))
    x.range <- max(processed$direct.shift) - min(processed$direct.shift)

    position <- position * x.range + min(processed$direct.shift)
    height <- height * y.range

    sfo1 <- get_parameter(nmrdata, 'sfo1', 'procs')
    width <- width * (x.range[2] - x.range[1]) * sfo1
  }

  #---------------------------------------
  # Defining a bound check function
  .check_bounds <- function(bounds) {

    if ( length(bounds) != 2 ) {
      err <- paste("All bounds must be vectors of two elements consisting",
                   "of a lower and upper bound.")
      stop(err)
    }

    if ( bounds[1] > bounds[2] ) {
      err <- paste("Lower bound must be smaller than upper bound.",
                   "Proceeding with current constraints will result in a",
                   "fit error.")
      warning(err)
    }

  }

  #---------------------------------------
  # Creating a list of bounds to loop through each in term
  bounds = list(position = position, height = height, width = width,
                fraction.gauss = fraction.gauss)

  for ( parameter in names(bounds) ) {
    if ( length(bounds[[parameter]]) > 0 ) {
      .check_bounds(bounds[[parameter]])
      lower[[parameter]] <- bounds[[parameter]][1]
      upper[[parameter]] <- bounds[[parameter]][2]
    }
  }

  # Fraction gauss is a little different because it must be 0-1
  lower$fraction.gauss[lower$fraction.gauss < 0] <- 0
  upper$fraction.gauss[upper$fraction.gauss > 0] <- 1

  # Ensuring that parameters are only widened if desired
  columns <- c('position', 'width', 'height', 'fraction.gauss')

  new.lower <- unlist(lower[ , columns])
  old.lower <- unlist(object@bounds$lower[ , columns])

  new.upper <- unlist(upper[ , columns])
  old.upper <- unlist(object@bounds$upper[ , columns])
  
  if (! widen ) {
    new.lower <- ifelse(new.lower < old.lower, old.lower, new.lower)
    new.upper <- ifelse(new.upper > old.upper, old.upper, new.upper)
  }

  object@bounds$lower[ , columns] <- new.lower
  object@bounds$upper[ , columns] <- new.upper

  validObject(object)
  object
})



#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRResonance1D",
  function(object, position = NULL, height = NULL, width = NULL, 
           relative = FALSE, widen = FALSE) {

  # Initializing bounds
  object <- .initialize_bounds(object)
  peaks <- object@peaks
  lower <- object@bounds$lower
  upper <- object@bounds$upper

  #---------------------------------------
  # Defining a bound check function
  .check_bounds <- function(bounds) {

    if ( length(bounds) != 2 ) {
      err <- paste("All bounds must be vectors of two elements consisting",
                   "of a lower and upper bound.")
      stop(err)
    }

    err2 <- "Proceeding with current constraints will result in a fit error."

    if ( bounds[1] > 0 ) {
      err <- paste("Lower offsets must be negative so that resulting bounds",
                   "include initial values.", err2)
      stop(err)
    }

    if ( bounds[2] < 0 ) {
      err <- paste("Upper offsets must be positive so that resulting bounds",
                   "include initial values.", err2)
      stop(err)
    }

    if ( bounds[1] > bounds[2] ) {
      err <- paste("Lower bound must be smaller than upper bound.", err2)
      warning(err)
    }

  }

  #---------------------------------------
  # Creating a list of bounds to loop through each in term
  bounds = list(position = position, height = height, width = width)

  for ( parameter in names(bounds) ) {
    if ( length(bounds[[parameter]]) > 0 ) {
      .check_bounds(bounds[[parameter]])

      lower.offset <- bounds[[parameter]][1]
      upper.offset <- bounds[[parameter]][2]
      
     if ( relative ) {
        lower.offset <- lower.offset*peaks[[parameter]]
        upper.offset <- upper.offset*peaks[[parameter]]
      } 

      lower[[parameter]] <- peaks[[parameter]] + lower.offset
      upper[[parameter]] <- peaks[[parameter]] + upper.offset
    }
  }

  # Ensuring that parameters are only widened if desired
  columns <- c('position', 'width', 'height', 'fraction.gauss')

  new.lower <- unlist(lower[ , columns])
  old.lower <- unlist(object@bounds$lower[ , columns])

  new.upper <- unlist(upper[ , columns])
  old.upper <- unlist(object@bounds$upper[ , columns])
  
  if (! widen ) {
    new.lower <- ifelse(new.lower < old.lower, old.lower, new.lower)
    new.upper <- ifelse(new.upper > old.upper, old.upper, new.upper)
  }

  object@bounds$lower[ , columns] <- new.lower
  object@bounds$upper[ , columns] <- new.upper

  validObject(object)
  object
})



#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRResonance1D",
  function(object, position = TRUE,  height = TRUE, width = TRUE, 
           nmrdata = NULL, widen = FALSE) { 

  # First, do a single pass over general bounds with no reference
  if ( height )  gen.height <- c(0, Inf)
  else gen.height <- NULL

  if ( width ) gen.width <- c(0.003, 3)
  else gen.width <- NULL

  object <- set_general_bounds(object, height = gen.height, width = gen.width,
                               widen = widen)

  # Adding position offsets
  if ( position ) {
    object <- set_offset_bounds(object, position = c(-0.1, 0.1), widen = widen)
  }

  # If nmrdata is provided, add further constraints  
  if (! is.null(nmrdata) ) {
    
    if ( class(nmrdata) != 'NMRData1D' ) {
      err <- '"nmrdata" must be a valid NMRData1D object.'
      stop(err)
    } else {
      validObject(nmrdata)
    }

    if ( position )  gen.position <- c(0, 1)
    else gen.position <- NULL

    if ( height ) gen.height <- c(0, 1.5)
    else gen.height <- NULL

    if ( width ) gen.width <- c(0, 0.2)
    else gen.width <- NULL

    object <- set_general_bounds(object, position = gen.position, 
                                 height = gen.height, width = gen.width,
                                 nmrdata = nmrdata, widen = widen)
  }

  object
  })



#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRResonance1D",
  function(object, peak.type) {

    # Initializing bounds
    object <- .initialize_bounds(object)
    peaks <- object@peaks
    lower <- object@bounds$lower
    upper <- object@bounds$upper

    # Getting rid of empty spaces and capitals
    peak.type <- tolower(gsub('\\s', '', peak.type))
    peak.types <- c('lorentz', 'voigt', 'gauss', 'any')
    peak.type <- pmatch(peak.type, peak.types)

    if ( peak.type == 1 ) {
      lower$fraction.gauss <- 0
      upper$fraction.gauss <- 0
      peaks$fraction.gauss <- 0
    } else if ( peak.type == 2 ) {
      lower$fraction.gauss <- 1e-6
      upper$fraction.gauss <- 1 - 1e-6

      logic <- peaks$fraction.gauss < lower$fraction.gauss
      peaks[logic , 'fraction.gauss'] <- lower$fraction.gauss[logic] + 1e-6
      logic <- peaks$fraction.gauss > upper$fraction.gauss
      peaks[logic , 'fraction.gauss'] <- upper$fraction.gauss[logic] - 1e-6
    } else if ( peak.type == 3 ) {
      lower$fraction.gauss <- 1
      upper$fraction.gauss <- 1
      peaks$fraction.gauss <- 1
    } else if ( peak.type == 4 ) {
      lower$fraction.gauss <- 0
      upper$fraction.gauss <- 1
    } else {
      peak.types <- paste(peak.types, collapse = ', ')
      err <- sprintf('Peak type must be one of %s', peak.types)
      stop(err)
    }

    object@bounds <- list(lower = lower, upper = upper)
    object@peaks <- peaks

    object
})

