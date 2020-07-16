# Definition of a super class structure for NMR data.



#------------------------------------------------------------------------------
#' Super class combining NMR data and scan parameters
#' 
#' Fourier transformed data is stored as a tibble (data.frame-like structure)
#' in the processed slot with the acquisition and processing parameters stored
#' as lists in separate slots. If the data is loaded from a JCAMP-DX file
#' without clear "acqus" and "procs" division, the parameters are stored in a
#' generic parameters list slot. Both parameter lists and chemical shift data
#' is clearly divided into "direct" or "indirect" components. If the data is
#' 1D, the "indirect" component is either omitted or NULL.
#' 
#' @slot processed A tibble containing chemical shift and intensity data.
#' @slot parameters A list of generic parameters, divided into "direct" and
#'                  "indirect" dimension components.
#' @slot acqus A list of parameters specific to acqus file, divided into
#'             "direct" and "indirect" dimension components.
#' @slot procs A list of parameters sepcific to procs file, divided into
#'             "direct" and "indirect" dimension components.
#' 
#' @name NMRData-class
#' @export
NMRData <- setClass("NMRData",
  representation(processed = "tbl",
                 parameters = "list",
                 acqus = "list",
                 procs = "list"),
  prototype(parameters = list(direct = list(), indirect = list()),
            acqus = list(direct = list(), indirect = list()),
            procs = list(direct = list(), indirect = list())))



#==============================================================================>
#  Validity check
#==============================================================================>



# A basic check of the components 
validNMRData <- function(object) {

  valid <- TRUE
  err <- c()

  processed <- object@processed
  acqus <- object@acqus
  procs <- object@procs

  # Dirty shortcut to avoid retyping the same thing
  add_err <- function(msg) {
    valid <- FALSE
    err <- c(err, msg)
  }

  #---------------------------------------
  # First the processed data

  # Processed should be a tibble
  msg <- '"processed" should be a tidyverse tibble'
  if (! "tbl" %in% class(processed) ) add_err(msg)

  # Intensity must be an ncmplx or ncmplx2 object
  msg <- '"intensity" columns must be "vctrs_cmplx1" or "vctrs_cmplx2" class'
  logic <- any(class(processed$intensity) %in% c('vctrs_cmplx1','vctrs_cmplx2'))
  if (! logic ) add_err(msg)

  if ( valid ) TRUE
  else err
}

setValidity("NMRData", validNMRData)



#==============================================================================>
#  Helper functions for 1D and 2D constructors
#==============================================================================>



#------------------------------------------------------------------------
#' Read Bruker acqus parameters
#' 
#' Reads the acqus parameters across all acquired dimensions. Data from each
#' dimension is stored as a separate list item, i.e., list(direct = [acqus],
#' indirect = [acqu2s]). Since the acqus files are in JCAMP format, the actual
#' parsing is just a thin wrapper around read_jcamp(), process_jcamp(), and
#' flatten_jcamp() that reads Bruker scan parameters and puts them in a flat
#' list.
#' 
#' @param path Character string that points to a scan directory.
#' @param ... Arguments passed into process_jcamp().
#' 
#' @return A list made up of nested lists with the processed acqus entries.
#' 
#' @export
read_acqus <- function(path, ...) {

  # Making sure that the path is a directory
  if (! dir.exists(path)) {
    msg <- '"path" must point to a scan directory containing acqus files.'
    stop(msg)
  }

  # Generating list of acqus file
  all.files <- list.files(path)

  # For now, valid files are restricted to 1D acqus and 2D acqu2s
  valid.files <- c('acqus'='direct', 'acqu2s'='indirect')
  acqus.files <- all.files[all.files %in% names(valid.files)]

  # Checking to make sure that at least some files were found
  err <- 'No acqus files found.'
  if (length(acqus.files) == 0) stop(err)

  # Combining with initial path and tacking on names
  acqus.paths <- file.path(path, acqus.files)
  names(acqus.paths) <- valid.files[acqus.files]

  # Defining read and process function
  f_process <- function(path, ...) {
    out <- read_jcamp(path, ...)
    flatten_jcamp(out)
  }

  # Generating nested list
  lapply(acqus.paths, f_process, ...)
}



#------------------------------------------------------------------------------
# Internal function for validating pdata path
.validate_pdata <- function(path, number) {

  err <- '"path" must point to an experiment directory containing pdata.'

  # First, check if current directory exists
  if (! dir.exists(path)) stop(err)

  # Directory must contain pdata
  logic <- ! 'pdata' %in% list.dirs(path, full.names = FALSE, recursive = FALSE)
  if ( logic ) stop(err)

  # pdata must contain folders
  pdata.path <- file.path(path, 'pdata')
  dirs <- list.dirs(pdata.path, full.names = FALSE, recursive = FALSE)
  
  err <- 'No directories found within pdata.'
  if ( length(dirs) == 0 ) stop(err)

  # Choosing default number if necessary
  if ( is.na(number) ) number <- dirs[1]
  
  file.path(path, 'pdata', number)

}



#------------------------------------------------------------------------
#' Read Bruker procs parameters
#' 
#' Reads the procs parameters across all acquired dimensions. Data from each
#' dimension is stored as a separate list item, i.e., list(direct = [procs],
#' indirect = [proc2s]). Since the procs files are in JCAMP format, the actual
#' parsing is just a thin wrapper around read_jcamp(), process_jcamp(), and
#' flatten_jcamp() that reads Bruker scan parameters and puts them in a flat
#' list.
#' 
#' @param path Character string that points to a scan directory.
#' @param number The processing file number to use. Defaults to smallest number
#'               in pdata folder.
#' @param dimension The dimension of the scan parameters (1 or 2). Defaults to
#'                  loading all
#' @param ... Arguments passed into process_jcamp().
#' 
#' @return A list made up of nested lists with the processed procs entries.
#' 
#' @export
read_procs <- function(path, number = NA, ...) {

  path <- .validate_pdata(path, number)

  # Generating list of procs file
  all.files <- list.files(path)

  # For now, valid files are restricted to 1D acqus and 2D acqu2s
  valid.files <- c('procs'='direct', 'proc2s'='indirect')
  procs.files <- all.files[all.files %in% names(valid.files)]

  # Checking to make sure that at least some files were found
  err <- 'No procs files found.'
  if (length(procs.files) == 0) stop(err)

  # Combining with initial path and tacking on names
  procs.paths <- file.path(path, procs.files)
  names(procs.paths) <- valid.files[procs.files]

  # Defining read and process function
  f_process <- function(path, ...) {
    out <- read_jcamp(path, ...)
    flatten_jcamp(out)
  }

  # Generating nested list
  lapply(procs.paths, f_process, ...)
}



#==============================================================================>
#  Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
#' Access specific parameter
#'
#' Looks for named parameter in one of either acqus or procs before defaulting
#' to the generic parameters list.
#'
#' @param NMRData object.
#' @param parameter Name of parameter.
#' @param location One of either "acqus" or "procs".
#' @param dimension One of either "direct" or "indirect".
#' @param error TRUE to issue an error if parameter not found.
#'
#' @name get_parameter
#' @export
setGeneric("get_parameter", 
  function(object, ...) standardGeneric("get_parameter"))

#' @rdname get_parameter
#' @export
setMethod("get_parameter", "NMRData", 
  function(object, parameter, location = "acqus", 
           dimension = "direct", error = FALSE) {

    out <- slot(object, location)[[dimension]][[parameter]]
    if ( is.null(out) ) out <- object@parameters[[parameter]]

    err <- sprintf('Parameter "%s" not found', parameter)
    if ( error && is.null(out) ) stop(err) 

    out
  
  })



#------------------------------------------------------------------------------
#' Access Fourier-transformed data
#' @name processed
#' @export
setGeneric("processed", 
  function(object, ...) standardGeneric("processed"))

#' @rdname processed
#' @export
setMethod("processed", "NMRData", 
  function(object) object@processed)

#' Set Fourier-transformed data
#' @name processed-set
#' @export
setGeneric("processed<-", 
  function(object, value) standardGeneric("processed<-"))

#' @rdname processed-set
#' @export
setReplaceMethod("processed", "NMRData",
  function(object, value) {
    object@processed <- value
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
#' Access generic parameters
#' @name parameters
#' @export
setGeneric("parameters", 
  function(object, ...) standardGeneric("parameters"))

#' @rdname parameters
#' @export
setMethod("parameters", "NMRData", 
  function(object) object@parameters)

#' Set generic parameters
#' @name parameters-set
#' @export
setGeneric("parameters<-", 
  function(object, value) standardGeneric("parameters<-"))

#' @rdname parameters-set
#' @export
setReplaceMethod("parameters", "NMRData",
  function(object, value) {
   object@parameters <- value
   validObject(object)
   object 
  })



#------------------------------------------------------------------------------
#' Access acqus parameters
#' @name acqus
#' @export
setGeneric("acqus", 
  function(object, ...) standardGeneric("acqus"))

#' @rdname acqus
#' @export
setMethod("acqus", "NMRData", 
  function(object) object@acqus)

#' Set acqus parameters
#' @name acqus-set
#' @export
setGeneric("acqus<-", 
  function(object, value) standardGeneric("acqus<-"))

#' @rdname acqus-set
#' @export
setReplaceMethod("acqus", "NMRData",
  function(object, value) {
   object@acqus <- value
   validObject(object)
   object 
  })



#------------------------------------------------------------------------------
#' Access procs parameters
#' @name procs
#' @export
setGeneric("procs", 
  function(object, ...) standardGeneric("procs"))

#' @rdname procs
#' @export
setMethod("procs", "NMRData", 
  function(object) object@procs)

#' Set procs parameters
#' @name procs-set
#' @export
setGeneric("procs<-", 
  function(object, value) standardGeneric("procs<-"))

#' @rdname procs-set
#' @export
setReplaceMethod("procs", "NMRData",
  function(object, value) {
   object@procs <- value
   validObject(object)
   object 
  })



#==============================================================================>
#  Values and processing
#==============================================================================>



#------------------------------------------------------------------------
#' Extract or calculate peak lineshape values
#' 
#' This is a core function used to convert an abstract NMRScaffold or NMRData
#' object into a tibble of chemical shifts and intensity values. The method for
#' this depends on the nature of the object. If the input is NMRData, the
#' output just extracts internal data for output. However, if the input is an
#' NMRScaffold object, new lineshape is data based on the chemical shifts and
#' sweep frequency provided. For 2D scaffolds, direct and indirect chemical
#' shifts will automatically generate a grid if unique.
#' 
#' @param object An NMRScaffold or NMRData object.
#' @param direct.shift Vector of chemical shift data in ppm.
#' @param indirect.shift Similar to direct.shift, but in the indirect dimension.
#' @param direct.sf Sweep frequency (in MHz) -- needed to convert peak widths
#'                  from Hz to ppm. In most cases, it is recommended to set a
#'                  single default value using nmroptions$direct$sf = ..., but
#'                  an override can be provided here.
#' @param indirect.sf Similar to direct.sf, but in the indirect dimension.
#' @param sum.level The level at which values should be calculated, with "all"
#'                  resulting in a single set of values from the sum of all
#'                  underlying components. Alternatively, separate values can be
#'                  generated at the "mixture", "species", and "resonance"
#'                  levels.
#' @param domain A character specifying which real/imaginary domains to output
#'               for intensity data. For 1D data "r/i" outputs both real and
#'               imaginary data, with "r" and "i" used to output one of either
#'               real or imaginary. For 2D data, this is expanded to
#'               "rr/ri/ir/ii". All domains are outputted by default.
#' @param use.cmplx1 cmplx1 is the internal representation of complex() that
#'                   matches cmplx2 for 2D data. It is replaced by default to
#'                   complex() for user-facing functions. However, this behavior
#'                   can be overriden here.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A tibble with "direct.shift", and "intensity" columns, as well as
#'         "indirect.shift" for 2D data and any combination of "resonance",
#'         "species", and "mixture" depending on the sum.level parameter and the
#'         input data.
#' 
#' @name values
#' @export
setGeneric("values", 
  function(object, ...) standardGeneric("values")
)



#------------------------------------------------------------------------
#' Add baseline to data
#' 
#' Baseline terms are defined in the same way as for NMRFit1D/NMRFit2D, so it
#' should be possible to apply these directly.
#' 
#' @param object An NMRData object.
#' @param baseline Vector of baseline spline parameters, where the domain of the
#'                 baseline must match the domain of the underlying data, i.e.,
#'                 a cmplx1 or complex vector is needed for 1D data and a cmplx2
#'                 vector is needed for 2D data. If a numeric vector is
#'                 provided, the same baseline is applied to all domains. Note
#'                 that for 2D data, the baseline components represent both
#'                 dimensions on an even grid and must therefore have a length
#'                 that is a perfect square. These values themselves represent
#'                 slices of baseline along the direct dimension in order of
#'                 increasing indirect dimension chemical shift.
#' @param knots A vector of relative internal knots ranging from 0
#'              (corresponding to lower direct.shift) and 1 (corresponding to
#'              upper direct.shift). Boundary knots are added automatically and
#'              should not be included.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A new NMRData object with modified intensity values.
#' 
#' @name add_baseline
#' @export
setGeneric("add_baseline", 
  function(object, ...) standardGeneric("add_baseline")
)



#------------------------------------------------------------------------
#' Add phase to data
#' 
#' Manual phase correction meant to apply phase terms calculated through
#' NMRFit. Note that phase terms are specified in radians unless "degrees" is
#' set to TRUE.
#' 
#' @param object An NMRData object.
#' @param phase A vector of phase correction terms -- 0 order, 1st order in the
#'              direct dimension, 1st order in the indirect.dimension.
#' @param degrees TRUE if phase terms are provided in degrees, FALSE by default.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A new NMRData object with modified intensity values.
#' 
#' @name add_phase
#' @export
setGeneric("add_phase", 
  function(object, ...) standardGeneric("add_phase")
)



#==============================================================================>
#  Defining list and data.frame like behaviour
#==============================================================================>



#------------------------------------------------------------------------------
#' Convert NMRData object to list
#' @export
as.list.NMRData <- function(x) {
  list(processed = x@processed, parameters = x@parameters, 
       procs = x@procs, acqus = x@acqus)
}

setMethod("as.list", "NMRData", as.list.NMRData)



#------------------------------------------------------------------------------
#' Convert NMRData object to a tidyverse tibble
#' @export
as_tibble.NMRData <- function(x) {
  x@processed
}



#==============================================================================>
#  Filter functions
#==============================================================================>



#------------------------------------------------------------------------------
#' Internal filter functions used for both direct and indirect shifts
.filter_shift <- function(x, lower, upper, round, align) {

  # If round is FALSE just do a basic filter and return
  logic <- (x > lower) & (x < upper)
  if ( round == FALSE ) return(logic)
  
  # Otherwise figure out rounding length and proceed from there
  n <- sum(logic)
  n.2 <- log(n)/log(2)
  if ( round == "up" ) n <- 2^ceiling(n.2)
  else if ( round == "down" ) n <- 2^floor(n.2)
  else {
    err <- '"round" must be one of "up", "down", or FALSE'
    stop(err)
  }

  # To ensure the following works, sort from smallest to largest
  y <- sort(x)

  # Align as required
  if ( align == "left" ) {
    logic <- y > lower
    start <- which(logic)[1]
    stop <- start + n - 1
  }
  else if ( align == "right" ) {
    logic <- y < upper
    stop <- which(logic)[sum(logic)]
    start <- stop - n + 1
  }
  else if ( align == "middle" ) {
    logic <- (y > lower) & (y < upper)

    n.interval <- sum(logic)
    start <- which(logic)[1]
    stop <- which(logic)[n.interval]

    n.mod <- n - n.interval
    start <- start - round(n.mod/2)
    stop <- stop + (n.mod - round(n.mod/2))
  } else {
    err <- '"align" must be one of "left", "right", or "middle"'
  }

  # Generating output
  msg <- "Not enough data points to maintain specified power of 2."
  if ( (start < 1) || (stop > length(y)) ) stop(msg)

  lower <- y[start]
  upper <- y[stop]

  (x > lower) & (x < upper)
}



#------------------------------------------------------------------------------
#' Filter NMRData based on chemical shift
#' 
#' Filter processed data to include only those points that are contained
#' between a set of lower and upper bounds on chemical shift. Data can be
#' filtered based on direct or indirect chemical shifts.
#' 
#' @param object An NMRData object.
#' @param lower A lower bound for chemical shift.
#' @param upper An upper bound for chemical shift.
#' @param round One of either "up", "down", or FALSE. Convolution is
#'              considerably faster when performed on vector lengths that are
#'              powers of two, and as such, it may be beneficial to round the
#'              total number remaining points, even if it means ignoring one or
#'              both of the bounds. Use this option to round the total number of
#'              filtered points up or down.
#' @param align One of either "lower", "upper", or "middle". Used to align
#'              points to the specified bounds, when using the "round" option
#'              prevents the specified bounds to be matched exaxctly.
#' 
#' @return An NMRData1D object with filtered processed data.
#' 
#' @name filter_direct
#' @export
setGeneric("filter_direct", 
  function(object, ...) {
    standardGeneric("filter_direct")
  })

#' @rdname filter_direct
#' @export
setMethod("filter_direct", "NMRData", 
  function(object, lower = NA, upper = NA, round = FALSE, align = "middle") {

    d <- object@processed 

    # Setting default bounds if not present
    if ( is.na(lower) ) lower <- min(d$direct.shift)
    if ( is.na(upper) ) upper <- max(d$direct.shift)

    # Directly apply internal filter function
    d <- filter(d, .filter_shift(direct.shift, lower, upper, round, align))
    object@processed <- d

    object
  })

#------------------------------------------------------------------------------
#' @rdname filter_direct
#' @export
setGeneric("filter_indirect", 
  function(object, ...) {
    standardGeneric("filter_indirect")
  })

#' @rdname filter_direct
#' @export
setMethod("filter_indirect", "NMRData", 
  function(object, lower = NA, upper = NA, round = FALSE, align = "middle") {

    d <- object@processed 

    if (! 'indirect.shift' %in% colnames(d) ) return(object)

    # Setting default bounds if not present
    if ( is.na(lower) ) lower <- min(d$indirect.shift)
    if ( is.na(upper) ) upper <- max(d$indirect.shift)

    # Apply internal filter function grouped by direct shift
    d <- group_by(d, direct.shift) %>%
           filter(d, .filter_shift(indirect.shift, lower, upper, round, align))
    object@processed <- ungroup(d)

    object
  })



#==============================================================================>
#  Processing
#==============================================================================>



#------------------------------------------------------------------------
#' Correct spectrum phase
#' 
#' This function applies phase correction to NMRData1D and NMRData2D objects.
#' The nature of the phase correction vector depends on the object. The 0 order
#' phase correction term is always taken as the phase correction at 0 ppm.
#' 
#' @param object An NMRData object.
#' @param phase A vector of phase corrections, from zero-order and up. 1st order
#'              correction for NMRData1D objects require a vector of length 2 (0
#'              order term, then first order term) while NMRData2D objects
#'              require a vector of length 3 (0 order term, first order in the
#'              direct dimension, 1st order in the indirect dimension).
#' @param degrees TRUE if the phase angles are in degrees, FALSE if in radians.
#
#' @return NMRData object with corrected phase.
#'
#' @name apply_phase
#' @export
setGeneric("apply_phase", 
  function(object, ...) standardGeneric("apply_phase")
)
