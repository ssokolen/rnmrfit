# Definition of a class structure for 1D NMR data.



#------------------------------------------------------------------------------
#' A class combining NMR data and scan parameters
#' 
#' An extension of the generic NMRData class to provide 1D-specific methods
#' such as constructors and plotting.
#' 
#' @rdname NMRData1D
#' @export
NMRData1D <- setClass("NMRData1D", contains = "NMRData")

# Validity testing consists of simply checking the processed data.frame columns
validNMRData1D <- function(object) {

  valid <- TRUE
  err <- c()

  processed <- object@processed
  acqus <- object@acqus
  procs <- object@procs

  # Checking processed columns
  valid.columns <- c('direct.shift', 'intensity')
  msg <- sprintf('Processed data may only have the following columns: %s',
                 paste(valid.columns, collapse = ', '))
  if (! identical(valid.columns, colnames(processed)) ) {
    valid <- FALSE
    err <- c(err, msg)
  }

  if ( valid ) TRUE
  else msg
}

setValidity("NMRData1D", validNMRData1D)



#==============================================================================>
#  Constructors and data loading
#==============================================================================>



#------------------------------------------------------------------------------
#' Constructors for generating an NMRData1D object
#' 
#' \code{nmrdata_1d()} can be used as a generic constructor method that handles
#' different types of input (currently limited to Bruker pdata directories and
#' JCAMP=DX files). If the input path points to a file, it is assumed to be a
#' JCAMP-DX file, and if the input path points to a directory, it is assumed to
#' be Bruker scan directory. The specific constructor functions can also be
#' called using \code{nmrdata_1d_from_pdata()} or
#' \code{nmrdata_1d_from_jcamp()}.
#' 
#' The ability to read data from a raw csv file is provided as a placeholder
#' until more data types are supported. csv files are assumed to have chemical
#' shift as 1st column, real domain data as 2nd column, and imaginary domain
#' data as 3rd column.
#' 
#' @param path Path to a Bruker scan directory or JCAMP-DX file.
#' @param procs.number Specifies pdata directory to load. Defaults to the lowest
#'                     available number. Ignored if loading from a JCAMP-DX
#'                     file.
#' @param blocks.number Specifies block number to use when loading from a JCAMP-
#'                      DX file. Defaults to the first block encountered.
#'                      Ignored if loading from a pdata directory.
#' @param ntuples.number Specifies ntuple entry number to use when loading from
#'                       a JCAMP-DX file. Defaults to the first ntuple entry
#'                       encountered. Ignored if loading from a pdata directory.
#' 
#' @return An NMRData1D object.
#' 
#' @export
nmrdata_1d <- function(path, procs.number = NA, 
                       blocks.number = 1, ntuples.number = 1) {

  # If path is a valid directory, treat as Bruker scan directory
  if ( file.exists(path) && dir.exists(path) ) {
    nmrdata_1d_from_pdata(path, procs.number) 
  }
  # If it's not a directory, it might be a file
  else if ( file.exists(path) )  {
    nmrdata_1d_from_jcamp(path, blocks.number, ntuples.number)
  }
  # Otherwise, error out
  else {
    err <- sprintf('Path "%s" does not point to a file or directory', path)
    stop(err)
  }

}



#------------------------------------------------------------------------------
#' @rdname nmrdata_1d
#' @export
nmrdata_1d_from_pdata <- function(path, procs.number = NA) {

  # First, loading procs parameters
  procs <- read_procs(path, procs.number)

  # Using the procs file to load the processed data
  processed <- read_processed_1d(path, procs$direct, procs.number)

  # Finally, loading the general acquisition parameters
  acqus <- read_acqus(path)

  # Returning class object
  new("NMRData1D", processed = processed, parameters = list(),
                   procs = procs, acqus = acqus)

}



#------------------------------------------------------------------------------
# Internal function for validating list parameters
.validate_param <- function(param.list, required.param) {

    missing <- ! required.param %in% names(param.list)

    err <- sprintf('The following parameters are missing: %s',
                    paste(required.param[missing], collapse=', '))
    if (any(missing)) stop(err)

    # Otherwise, returning the required parameters
    param.list[required.param]

}



#------------------------------------------------------------------------
#' Read 1D Bruker 1r/1i files
#' 
#' Reads processed bruker 1D files based on specified parameters.
#' 
#' @param path Character string that points to a Bruker scan directory.
#' @param procs.list A list of procs parameters that contains 'sw.p', 'si',
#'                   'sf', 'reverse', and 'offset' entries. This list can be
#'                   generated using read_procs() or through other means.
#' @param number The processing file number. Defaults to smallest number in
#'               pdata directory.
#' 
#' @return A tidyverse data_frame made of two columns -- "direct.shift" and
#'         "intensity", corresponding to direct dimension chemical shift and the
#'         complex spectrum intensity data stored in a cmplx1 vector,
#'         respectively.
#' 
#' @export
read_processed_1d <- function(path, procs.list, number = NA) {

  # Checking for required entries
  required.procs <- c('sw.p', 'si', 'sf', 'reverse', 'offset')
  procs.list <- .validate_param(procs.list, required.procs)

  # Extracting parameters
  sw.p <- procs.list$sw.p
  si <- procs.list$si
  sf <- procs.list$sf
  rv <- procs.list$reverse 
  ofs <- procs.list$offset

  # Doing some basic validation
  path <- .validate_pdata(path, number)
  
  # Setting file path 
  path.real <- file.path(path, '1r')
  path.imag <- file.path(path, '1i')

  # Reading binary files
  real.data <- safe_read(path.real, 'bin', size = 4, what = 'integer', n = si)
  imag.data <- safe_read(path.imag, 'bin', size = 4, what = 'integer', n = si)

	if ( (length(real.data) < si) | (length(imag.data) < si)){
    msg <- sprintf('Error reading 1r/1i files, file size does not match data')
		stop(msg)
	}

  # Adding a sign change to the imaginary domain to account for flip when
  # compared to raw fft of signal -- not sure why this is necessary...
  intensity  <- cmplx1(r = real.data, i = -imag.data)

  # Formatting the x-axis
  direct.shift <- seq(ofs, ofs - sw.p/sf, length.out = si)
  if (rv == 'yes')  direct.shift <- rev(direct.shift)

  # Combining output
  tibble(direct.shift = direct.shift, intensity = intensity)
}



#------------------------------------------------------------------------------
#' @rdname nmrdata_1d
#' @export
nmrdata_1d_from_rs2d <- function(path, procs.number = NA) {

  err <- '"path" must point to an experiment directory containing Proc'

  # Directory must contain Proc
  logic <- ! 'Proc' %in% list.dirs(path, full.names = FALSE, recursive = FALSE)
  if ( logic ) stop(err)

  # pdata must contain folders
  procs.path <- file.path(path, 'Proc')
  dirs <- list.dirs(procs.path, full.names = FALSE, recursive = FALSE)
  
  err <- 'No directories found within Proc.'
  if ( length(dirs) == 0 ) stop(err)

  # Choosing default number if necessary
  if ( is.na(procs.number) ) procs.number <- dirs[1]
  
  procs.path <- file.path(path, 'Proc', procs.number)

  # Picking off required parameters
  header.path <- file.path(procs.path, "header.xml")
  d <- xml2::read_xml(header.path)

  required.procs <- c(
    "SPECTRAL_WIDTH", 
    "MATRIX_DIMENSION_1D", 
    "OBSERVED_FREQUENCY", 
    "OFFSET_FREQ_1"
  )

  string.xpath <-"/header/params/entry/key[text()='%s']/../value/value"
  value.xpath <- sprintf(string.xpath, required.procs)
  procs.values <- vapply(value.xpath,  
    function(x) xml_text(xml_find_all(d, x)), ""
  )

  # Extracting parameters
  sw.p <- as.numeric(procs.values[1])
  si <- as.numeric(procs.values[2])
  sf <- as.numeric(procs.values[3])/1e6
  ofs <- as.numeric(procs.values[4])

  procs.list <- list(sw.p = sw.p, si = si, sf = sf, ofs = ofs)

  # Formatting the x-axis
  direct.shift <- seq(-sw.p/2/sf, sw.p/2/sf, length.out = si) - ofs/sf

  # Reading binary files
  binary.path <- file.path(procs.path, "data.dat")
  all.data <- safe_read(binary.path, 'bin', size = 4, 
                        what = 'double', n = si*2, endian = "big")
  
  real.data <- all.data[seq(1, si*2, by = 2)]
  imag.data <- all.data[seq(2, si*2, by = 2)]

  intensity <- cmplx1(r = real.data, i = imag.data)

  # Finally, combine the data
  processed <- tibble(direct.shift = direct.shift, intensity = intensity)

  # Returning class object
  new("NMRData1D", processed = processed, parameters = procs.list,
                   procs = list(), acqus = list()) 
}

#------------------------------------------------------------------------------
#' @rdname nmrdata_1d
#' @export
nmrdata_1d_from_jcamp <- function(path, blocks.number = 1, ntuples.number = 1) {

  jcamp <- read_jcamp(path, process.tags = TRUE, process.entries = TRUE)

  # Double check that specified block exists
  err <- sprintf("Specified block number not found in %s", path)
  if ( blocks.number > length(jcamp$blocks) ) stop(err)

  jcamp$blocks <- jcamp$blocks[blocks.number]

  # Check to make sure that ntuples exist
  err <- 'Import from JCAMP file currently limited to NTUPLES entries'
  if (! 'ntuples' %in% names(jcamp$blocks[[1]]) ) stop(err)

  # Double check that specified ntuple exists
  err <- sprintf('Specified NTUPLES number not found in block %i of %s', 
                 ntuples.number, path)
  if ( ntuples.number > length(jcamp$blocks[[1]]) ) stop(err)

  jcamp$blocks[[1]]$ntuples <- jcamp$blocks[[1]]$ntuples[ntuples.number]

  # Flattening file
  jcamp.flat <- flatten_jcamp(jcamp)

  # Extracting processed data from ntuples
  descriptors <- jcamp.flat$ntuples$descriptors
  pages <- jcamp.flat$ntuples$pages

  # Checking variable names
  variables <- tolower(descriptors$var.name)
  real.index <- which(str_detect(variables, 'spectrum.*real'))
  imag.index <- which(str_detect(variables, 'spectrum.*imag'))

  # Double check that the first entry is frequency
  err <- 'Import from JCAMP file currently limited to frequency abscissa'
  if (! str_detect(variables[1], 'freq') ) stop(err)

  # If both real and imaginary data isn't there, abort
  err <- 'Import from JCAMP file currently limited to real/imaginary spectra'
  if ( length(c(real.index, imag.index)) < 2 ) stop(err)

  # Proceed to extract data
  real.data <- pages[[real.index - 1]]
  imag.data <- pages[[imag.index - 1]]

  # Checking that frequency is the same
  real.frequency <- real.data[, 1]
  imag.frequency <- imag.data[, 1]
  err <- 'Mismatch in real and imaginary frequency, likely parsing error'
  if ( any(abs(real.frequency - imag.frequency) > 1e-6) ) stop(err)

  # Starting with raw values
  frequency <- real.frequency
  real.data <- real.data[, 2]
  imag.data <- imag.data[, 2]

  # Scaling if factors are non zero
  scale <- descriptors$factor[1]
  if (scale > 1e-6) frequency <- frequency*scale

  scale <- descriptors$factor[real.index]
  if (scale > 1e-6) real.data <- real.data*scale

  scale <- descriptors$factor[imag.index]
  if (scale > 1e-6) imag.data <- imag.data*scale

  # Offsetting if max-min difference is non zero
  d.max <- descriptors$max[1]
  d.min <- descriptors$min[1]
  if ( (d.max - d.min) > 1e-6 ) {
    frequency <- frequency - max(frequency) + d.max
  }

  d.max <- descriptors$max[real.index]
  d.min <- descriptors$min[real.index]
  if ( (d.max - d.min) > 1e-6 ) {
    real.data <- real.data - max(real.data) + d.max
  }

  d.max <- descriptors$max[imag.index]
  d.min <- descriptors$min[imag.index]
  if ( (d.max - d.min) > 1e-6 ) {
    imag.data <- imag.data - max(imag.data) + d.max
  }

  # Doing one final check on the direct shift to check on offset
  direct.shift <- frequency/jcamp.flat$sf
  
  delta <- jcamp.flat$offset - max(direct.shift)
  if ( abs(delta) > 1e-6 ) direct.shift <- direct.shift + delta

  # Finally, combine the data
  intensity <- cmplx1(r = real.data, i = imag.data)
  processed <- tibble(direct.shift = direct.shift, intensity = intensity)

  # Returning class object
  new("NMRData1D", processed = processed, parameters = jcamp.flat,
                   procs = list(), acqus = list()) 
}



#------------------------------------------------------------------------
#' @rdname nmrdata_1d
#' @export
nmrdata_1d_from_csv <- function(path) {

  d <- read_csv(path, col_names = FALSE, col_types = 'ddd')
  direct.shift <- pull(d, 1)
  intensity <- cmplx1(r = pull(d, 2), i = pull(d, 3))

  processed <- tibble(direct.shift = direct.shift, intensity = intensity)

  new("NMRData1D", processed = processed) 
}



#------------------------------------------------------------------------------
#' Convert NMRScaffold1D to NMRData1D
#' 
#' Generate lineshape from specified peaks.
#' 
#' @param object An NMRScaffold1D object
#' @param direct.shift Vector of chemical shift data in ppm.
#' 
#' @return An NMRData1D object.
#' 
#' @export
nmrdata_1d_from_scaffold <- function(object, direct.shift = NULL) {

  # Using the procs file to load the processed data
  processed <- values(object, direct.shift = direct.shift, use.cmplx1 = TRUE)

  # acqus just contains the sf
  acqus <- list(direct = list(sfo1 = object@sf))

  # Returning class object
  new("NMRData1D", processed = processed, acqus = acqus)
}



#==============================================================================>
#  Corrections
#==============================================================================>



#------------------------------------------------------------------------------
#' Correct phase
#' 
#' Applies basic entropy minimization algorithm to correct overall phase.
#' 
#' @param object An NMRData object.
#' @param iterations How many times to iterate the correction. A scaling term
#'                   must be calculated at the beginning of the correction and
#'                   it makes sense to update this scaling term once a rough
#'                   correction has been established -- so a minimum of 2
#'                   iterations is recommended.
#' 
#' @return A new NMRData object with corrected phase. If return.phase is TRUE, a
#'         two element vector is returned representing the 0 and 1st order phase
#'         terms.
#' 
#' @name correct_phase
#' @export
setGeneric("correct_phase", 
  function(object, ...) standardGeneric("correct_phase")
)



#' @rdname correct_phase
#' @export
setMethod("correct_phase", "NMRData1D", 
  function(object, iterations = 2) {

    # Choosing scaling factor
    d <- values(object, "r")
    x <- d$direct.shift
    y <- d$intensity 

    y <- y - median(y)

    r <- sgolayfilt(y, p = 3, n = 5, m = 2)
    h <- abs(r)/max(abs(r))
    sw <- max(x) - min(x)

    entropy <- sum(h*log(h))
    penalty <- sum(y[y < 0]^2)
    scaling <- abs(entropy/penalty)
    
    f <- function(p) {
      corrected <- add_phase(object, p)
      
      y <- values(corrected, "r")$intensity
      y <- y - median(y)
      r <- sgolayfilt(y, p = 3, n = 5, m = 2)
      h <- abs(r)/max(abs(r))

      entropy <- sum(h*log(h))
      penalty <- scaling*sum(y[y < 0]^2)

      -entropy + penalty
    }

    opt <- optim(c(0, 0), f, method = "L-BFGS-B", 
                 lower = c(-pi, -pi/sw), upper = c(pi, pi/sw))

    # Iterate if required
    out <- add_phase(object, opt$par) 

    if ( iterations > 1 ) {
      out <- correct_phase(out, iterations = iterations - 1)
    }

    out
  })



#------------------------------------------------------------------------------
#' Correct reference
#' 
#' A very basic algorithm used to correct the chemical shift scale. The biggest
#' positive peak within a certain distance of the the current 0 ppm region is
#' set as the new 0.
#' 
#' @param object An NMRData object.
#' @param span The ppm distance to search for the biggest peak starting from the
#'             current value of 0 ppm.
#' 
#' @return A new NMRData object with corrected chemical shift data.
#' 
#' @name correct_reference
#' @export
setGeneric("correct_reference", 
  function(object, ...) standardGeneric("correct_reference")
)



#' @rdname correct_reference
#' @export
setMethod("correct_reference", "NMRData1D", 
  function(object, span = 1) {

    # Choosing scaling factor
    d <- values(object, "r")
    x <- d$direct.shift
    y <- d$intensity

    logic <- (x > -span) & (x < span)
    x <- x[logic]
    y <- y[logic]

    x.shift <- x[y == max(y)]

    d <- processed(object)
    d$direct.shift <- d$direct.shift - x.shift[1]
    object@processed <- d

    object
  })




#==============================================================================>
#  Processing
#==============================================================================>



#------------------------------------------------------------------------
# Extracting lineshape values
#' @rdname values
#' @export
setMethod("values", "NMRData1D", 
  function(object, domain = "r/i", use.cmplx1 = FALSE) {
    
    d <- object@processed
    d$intensity <- domain(d$intensity, domain, use.cmplx1)

    d
  })



#------------------------------------------------------------------------
# Add baseline
#' @rdname add_baseline
#' @export
setMethod("add_baseline", "NMRData1D", 
  function(object, baseline, knots = c()) {

    # Ensuring there are enough parameters
    err <- '"baseline" must be 1 or more elements longer than "knots"'
    if ( length(baseline) <= length(knots) ) stop(err)

    # Extracting processed values
    d <- values(object)

    n <- nrow(d)
    n.baseline <- length(baseline)

    # Generating baseline parameters based on input
    if ( "numeric" %in% class(baseline) ) {
      p <- c(baseline, baseline)
    } else {
      p <- c(Re(baseline), Im(baseline))
    }

    # Tacking out bounds to internal knots
    knots <- sort(c(knots, c(0, 1)))

    # Direct shift has to be scaled since knots are relative
    x.range <- range(d$direct.shift)
    x.span <- x.range[2] - x.range[1]

    # Feeding into Rust wrapper
    out <- .Call(
      "baseline_1d_wrapper", 
      x = as.double((d$direct.shift - x.range[1])/x.span), 
      y = as.double(rep(0, n*2)), 
      knots = as.double(knots),
      p = as.double(p), 
      n = as.integer(n),
      nb = as.integer(n.baseline),
      nk = as.integer(length(knots))
    )

    index.re <- 1:n
    index.im <- index.re + n

    d$intensity <- d$intensity + cmplx1(r = out[index.re], i = out[index.im])
    object@processed <- d

    object
  })



#------------------------------------------------------------------------
# Add phase
#' @rdname add_phase
#' @export
setMethod("add_phase", "NMRData1D", 
  function(object, phase, degrees = FALSE) {

    if ( length(phase) == 0 ) return(object)

    # Ensuring correct length
    err <- '"phase" must be of length 1 (0 order) or 2 (1st order)'
    if (! length(phase) %in% c(1, 2) ) stop(err)

    # Extracting processed values
    d <- values(object)

    n <- nrow(d)
    n.phase <- length(phase)

    # Modifying phase angles if required
    if ( degrees ) phase <- phase/180*pi

    # Feeding into Rust wrapper
    out <- .Call(
      "phase_1d_wrapper", 
      x = as.double(d$direct.shift), 
      y = as.double(c(Re(d$intensity), Im(d$intensity))), 
      p = as.double(phase), 
      n = as.integer(n),
      np = as.integer(n.phase)
    )

    index.re <- 1:n
    index.im <- index.re + n

    d$intensity <- cmplx1(r = out[index.re], i = out[index.im])
    object@processed <- d

    object
  })



#------------------------------------------------------------------------
# Add noise
#' @rdname add_noise
#' @export
setMethod("add_noise", "NMRData1D", 
  function(object, sd = 0.01) {

    # Extracting processed values
    d <- values(object)
    
    sd <- sd*max(Re(d$intensity))
    n <- nrow(d)

    d$intensity <- d$intensity + 
      cmplx1(r = rnorm(n, sd = sd), i = rnorm(n, sd = sd))
    object@processed <- d

    object
  })



#==============================================================================>
#  Formatting and printing
#==============================================================================>



#' @export
format.NMRData1D <- function(x, ...) {
  d <- processed(x)
  components <- paste(colnames(as_tibble(d$intensity)), collapse = ', ') 
  shift.range <- range(d$direct.shift)
  shift.range <- sprintf("%.3f ppm to %.3f ppm", shift.range[1], shift.range[2])
  sprintf('NMRData1D object (%s), %s\n', components, shift.range)
}

#' @export
print.NMRData1D <- function(x, ...) cat(format(x))

#' @export
setMethod("show", "NMRData1D", 
  function(object) cat(format(object))
  )

#' @export
summary.NMRData1D <- function(object, ...) {
  d <- object@processed
  summary(bind_cols(direct.shift = d$direct.shift, 
                    as_tibble(d$intensity)))
}

#' @export
range.NMRData1D <- function(object, ...) {
  d <- object@processed
  list(direct.shift = range(d$direct.shift),
       intensity = range(Re(d$intensity)))
}

#' @export
is_vector_s3.NMRData1D <- function(x) FALSE

#' @export
type_sum.NMRData1D <- function(x) "NMRData1D"

#' @export
obj_sum.NMRData1D <- function(x) format(x)



#==============================================================================>
#  Plotting
#==============================================================================>

#------------------------------------------------------------------------------
#' Plot NMRData1D object
#' 
#' Convenience function that generates a plot of the spectral data.
#' 
#' @param x An NMRData1D object.
#' @param domain One of either 'r' or 'i' corresponding to either real or
#'               imaginary data. 
#' @param legendgroup Unique value for grouping legend entries.
#' @param color Intended for internal use -- line colour as character 
#' @param name Intended for internal use -- data name as character.
#' 
#' @return A plot_ly plot.
#' 
#' @export
plot.NMRData1D <- function(x, domain = 'r', legendgroup = 1,
                           color = NULL, name = NULL) {

  legend.opts <- list(orientation = 'h', xanchor = "center", x = 0.5)

  f <- list(
    family = "Courier New, monospace",
    size = 16,
    color = "#7f7f7f"
  )

  xaxis <- list(
    title = "Direct chemical shift (ppm)",
    titlefont = f,
    autorange = "reversed"
  )
  yaxis <- list(
    title = "Intensity",
    titlefont = f
  )

  #---------------------------------------

  err <- '"domain" must be one of "r" or "i"'
  if (! domain %in% c("r", "i") ) stop(err)

  d <- values(x, domain = domain)
  x <- d$direct.shift
  y <- d$intensity

  if ( is.null(color) ) color <- "black"
  if ( is.null(name) ) name <- "Raw data"

  p <- plot_ly(x = x, y = y, color = I(color), 
               name = I(name), type = 'scatter', mode = 'lines',
               legendgroup = legendgroup) %>%
       layout(legend = legend.opts,
              xaxis = xaxis, yaxis = yaxis)

  p
}

setMethod("plot", "NMRData1D", plot.NMRData1D)



#------------------------------------------------------------------------------
#' Add to existing NMRData1D plot
#' 
#' Meant to be primarily an internal function that adds a new line to an
#' existing plot.
#' 
#' @param x An NMRData1D object.
#' @param p A plot_ly object.
#' @param domain One of either 'r' or 'i' corresponding to either real or
#'               imaginary data. are displayed in separate subplots.
#' @param legendgroup Unique value for grouping legend entries.
#' @param color Line colour as character
#' @param name Data name as character.
#' 
#' @return A plot_ly plot.
#' 
#' @export
lines.NMRData1D <- function(x, p, domain = 'r', legendgroup = 2,
                            color = NULL, name = NULL) {

  stop <- '"domain" must be one of "r" or "i"'
  if (! domain %in% c("r", "i") ) stop(err)

  d <- values(x, domain = domain)
  x <- d$direct.shift
  y <- d$intensity

  if ( is.null(color) ) color <- "black"
  if ( is.null(name) ) name <- "Raw data"

  p <- p %>%
    add_trace(x = x, y = y, color = I(color), 
              name = I(name), type = 'scatter', mode = 'lines',
              legendgroup = legendgroup)

  p
}

setMethod("lines", "NMRData1D", lines.NMRData1D)
