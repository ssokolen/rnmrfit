# Definition of a class structure for 1D resonance data.



#==============================================================================>
#  NMRResonance1D -- peak description 
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR resonance.
#' 
#' Essentially, this class is used to define coupling relationships to group
#' individual peaks into resonances. A resonance is a series of peaks that are
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
    name = 'character',
    id = 'character',
    sf = 'numeric',
    peaks = 'data.frame',
    couplings = 'data.frame',
    couplings.leeway = 'list',
    lower.bounds = 'data.frame',
    upper.bounds = 'data.frame'
  ),
  prototype = prototype(
    name = 'resonance',
    id = 'resonance',
    sf = numeric(0),
    couplings = data.frame(),
    couplings.leeway = list(position = 0, width = 0, 
                            fraction.gauss = 0, area = 0),
    lower.bounds = data.frame(),
    upper.bounds = data.frame()
  )
)



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
#' @param constant The coupling constant of the split (in ppm).
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
#' @param position TRUE to enforce position differences.
#' @param width TRUE to enforce constant width. 
#' @param fraction TRUE to enforce constant fraction gauss. 
#' @param area TRUE to enforce area ratios. 
#' 
#' @return An NMRResonance1D object with a new set of coupling constraints.
#' 
#' @export
enforce_couplings_1d <- function(nmrresonance, position = TRUE, width = TRUE,
                                 fraction = TRUE, area = TRUE) {

  # Checking validity
  if ( class(nmrresonance) != 'NMRResonance1D' ) {
    err <- '"nmrresonance" must be a valid NMRResonance1D object.'
    stop(err)
  }

  # Filtering peaks
  d <- nmrresonance@peaks

  # There must be more than one peak to add couplings
  if ( nrow(d) <= 1 ) return(nmrresonance)

  # Calculating
  index.1 <- 1:(nrow(d) - 1)
  index.2 <- 2:nrow(d)

  # Initializing NA values
  differences <- NA
  ratios <- NA

  if ( position ) differences <- d$position[index.2] - d$position[index.1]
  if ( area ) ratios <- d$height[index.2]/d$height[index.1]
  
  if (! width ) width <- NA
  if (! fraction ) fraction <- NA

  couplings <- data.frame(peak.1 = d$peak[index.1], peak.2 = d$peak[index.2],
                          position.difference = differences,
                          area.ratio = ratios, width = width, fraction = fraction)

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
#' @param sf Sweep frequency (in MHz) -- needed to convert coupling
#'                  constants from Hz to ppm. In most cases, it is recommended
#'                  to set a single default value using nmroptions$direct$sf =
#'                  ..., but an override can be provided here.
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
#'                              would enable one peak to be 24\% gauss while
#'                              another peak within the same resonance is 34\%.
#' @param area.leeway Similar to position.leeway but for peak areas. Determines
#'                    how strictly the coupling area ratios are enforced.
#' 
#' @return An NMRResonance1D object.
#' 
#' @export
nmrresonance_1d <- function(peaks, sf = nmroptions$direct$sf, id = NULL, 
                            width = 1, fraction.gauss = 0, 
                            position.leeway = 0, width.leeway = 0, 
                            fraction.gauss.leeway = 0, area.leeway = 0) {

  # Checking to make sure that sweep frequency is defined
  err <- '"sf" must be provided as input or set nmroptions$direct$sf'
  if ( is.null(sf) ) stop(err)

  #---------------------------------------
  # Building peak list

  # Couplings are not added in every case
  add.couplings <- FALSE

  # If peaks is just a single number, it's better processed as a singlet
  if ( is.numeric(peaks) && ( length(peaks) == 1) ) peaks <- paste(peaks, 's')

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
  else if ( is.numeric(peaks) ) {
    add.constraints <- FALSE
    middle <- (max(peaks) + min(peaks))/2 
    range <- paste(min(peaks), '..', max(peaks), sep = '')
    if ( is.null(id) ) id <- paste(middle, 'm', range)
    peaks <- data.frame(peak = 1:length(peaks), 
                        position = peaks, width = width, height = 1, 
                        fraction.gauss = fraction.gauss)
  } else {
    err <- '"peaks" are defined from strings or numeric singlet ppm values.'
    stop(err)
  }

  #---------------------------------------
  # Adding coupling definitions

  # Starting with blanks first
  couplings <- data.frame()
  couplings.leeway = list(position = position.leeway, width = width.leeway,
                          fraction.gauss.leeway = fraction.gauss.leeway, 
                          area = area.leeway)

  #---------------------------------------
  # Adding infinite bounds

  # Handling bounds if they exist
  bounds <- list(lower = NULL, upper = NULL)

  # Selecting default values
  values <- list(lower = -Inf, upper = +Inf)
  columns <- c('position', 'width', 'height', 'fraction.gauss')

  for ( name in names(bounds) ) {
    peaks.copy <- peaks
    peaks.copy[ , columns] <- values[[name]]

    bounds[[name]] <- peaks.copy
  }

  # Fraction gauss must stay within 0-1
  bounds$lower$fraction.gauss <- 0
  bounds$upper$fraction.gauss <- 1

  nmrresonance = new('NMRResonance1D', id = id, peaks = peaks, sf = sf,
                                       couplings = couplings, 
                                       couplings.leeway = couplings.leeway,
                                       lower.bounds = bounds$lower,
                                       upper.bounds = bounds$upper)

  # And then updating if necessary
  if ( add.couplings ) enforce_couplings_1d(nmrresonance)
  else nmrresonance
}
