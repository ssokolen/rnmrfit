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



#==============================================================================>
#  Processing
#==============================================================================>



#------------------------------------------------------------------------
# Phase

#' @rdname apply_phase
#' @export
setMethod("apply_phase", "NMRData1D", 
  function(object, phase, degrees = TRUE) {

  # Summing up phase components
  n.phase <- length(phase)
  if ( n.phase == 0 ) return(object)

  # Changing phase angle units if required
  if (degrees) phase <- phase*pi/180

  n <- length(spectrum)
  phase.total <- rep(phase[1], n)

  processed <- processed(object)

  if ( length(phase) == 1) {
    phase.total <- phase
  } else if ( length(phase) == 2 ) {
    direct.shift <- processed$direct.shift
    phase.total <- phase.total + phase[2]*direct.shift
  } else {
    err <- '"phase" must be a vector of length 1 or 2.'
    stop(err)
  }

  im <- Im(processed$intensity)
  re <- Re(processed$intensity)

  intensity <- cmplx1(r = re * cos(phase.total) + im * sin(phase.total),
                      i= im * cos(phase.total) - re * sin(phase.total))

  processed$intensity <- intensity
  object@processed <- processed
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
#' @param components One of either 'r', 'i', or 'r/i' to include real, imaginary
#'                   or both components. If both components are selected, they
#'                   are displayed in separate subplots.
#' 
#' @return A plot_ly plot.
#' 
#' @export
plot.NMRData1D <- function(x, components = 'r') {

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
  d <- x@processed
  direct.shift <- d$direct.shift
  y.data <- d$intensity

  # Defining generic plot function
  f_init <- function(y, color, name) {
    p <- plot_ly(x = direct.shift, y = y, color = I(color), 
                 name = I(name), type = 'scatter', mode = 'lines',
                 legendgroup = 1) %>%
         layout(legend = legend.opts,
                xaxis = xaxis, yaxis = yaxis)
  }

  # Initializing the plot list
  plots <- list()

  # Checking which components to plot
  re <- grepl('r', components)
  im <- grepl('i', components)

  # Plotting 
  if ( re ) plots$r <- f_init(Re(y.data), 'black', 'Real')
  if ( im ) plots$i <- f_init(Im(y.data), 'grey', 'Imaginary')

  if ( length(plots) == 0 ) NULL
  else subplot(plots, shareX = TRUE, shareY = TRUE, 
               nrows = min(length(plots), 2))
}

setMethod("plot", "NMRData1D", plot.NMRData1D)
