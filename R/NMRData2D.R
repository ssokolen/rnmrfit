# Definition of a class structure for 2D NMR data.



#------------------------------------------------------------------------------
#' A class combining NMR data and scan parameters
#' 
#' An extension of the generic NMRData class to provide 2D-specific methods
#' such as constructors and plotting.
#' 
#' @rdname NMRData2D
#' @export
NMRData1D <- setClass("NMRData2D", contains = "NMRData")

# Validity testing consists of simply checking the processed data.frame columns
validNMRData2D <- function(object) {

  valid <- TRUE
  err <- c()

  processed <- object@processed
  acqus <- object@acqus
  procs <- object@procs

  # Checking processed columns
  valid.columns <- c('direct.shift', 'indirect.shift', 'intensity')
  msg <- sprintf('Processed data may only have the following columns: %s',
                 paste(valid.columns, collapse = ', '))
  if (! identical(valid.columns, colnames(processed)) ) {
    valid <- FALSE
    err <- c(err, msg)
  }

  if ( valid ) TRUE
  else msg
}

setValidity("NMRData2D", validNMRData2D)



#==============================================================================>
#  Constructors and data loading
#==============================================================================>



#------------------------------------------------------------------------------
#' Constructors for generating an NMRData1D object
#' 
#' \code{nmrdata_2d()} can be used as a generic constructor method that will
#' eventually handle different types of input (currently limited to Bruker
#' pdata directories). The specific constructor functions can also be called
#' using \code{nmrdata_2d_from_pdata()}.
#' 
#' @param path Path to a Bruker scan directory.
#' @param procs.number Specifies pdata directory to load. Defaults to the lowest
#'                     available number.
#' 
#' @return An NMRData1D object.
#' 
#' @export
nmrdata_2d <- function(path, procs.number = NA) {

  # If path is a valid directory, treat as Bruker scan directory
  if ( file.exists(path) && dir.exists(path) ) {
    nmrdata_2d_from_pdata(path, procs.number) 
  }
  # Otherwise, error out
  else {
    err <- sprintf('Path "%s" does not point to a file or directory', path)
    stop(err)
  }

}



#------------------------------------------------------------------------------
#' @rdname nmrdata_2d
#' @export
nmrdata_2d_from_pdata <- function(path, procs.number = NA) {

  # First, loading procs parameters
  procs <- read_procs(path, procs.number)

  # Using the procs file to load the processed data
  processed <- read_processed_2d(path, procs, procs.number)

  # Finally, loading the general acquisition parameters
  acqus <- read_acqus(path)

  # Returning class object
  new("NMRData2D", processed = processed, parameters = list(),
                   procs = procs, acqus = acqus)

}



#------------------------------------------------------------------------
#' Read 2D Bruker rr/ri/ir/ii files
#' 
#' Reads processed bruker 2D files based on specified parameters. As the
#' processed data can vary considerably based on the quadrature method used,
#' the resulting output will also vary a bit. If all 2rr/2ri/2ir/2ii files are
#' present, they are stored as a cmplx2 object that encodes all 4 domains. If
#' only two files are present, they are stored as a cmplx1 object with familiar
#' r/i domains. If no imaginary data is provided at all, intensity is left as a
#' cmplx1 object with an imaginary component set to zero.
#' 
#' @param path Character string that points to a Bruker scan directory.
#' @param procs.list A list of lists containing procs parameters with 'sw.p',
#'                   'si', 'sf', 'reverse', and 'offset'  entries for each of
#'                   the 'direct' and 'indirect' sublists. This list can be
#'                   generated using read_procs().
#' @param number The processing file number. Defaults to the smallest number in
#'               the pdata directory.
#' 
#' @return A tidyverse data_frame made up of three columns -- "direct.shift"
#'         containing the direct dimension chemical shift, "indirect.shift"
#'         containing the indirect dimension chemical shift, and "intensity"
#'         containing either a cmplx1 or cmplx2 vector.
#' 
#' @export
read_processed_2d <- function(path, procs.list, number = NA) {

  # The procs.list must contain appropriate sublists 
  err <- 'procs.list must contain two sublists named "direct" and "indirect".'
  logic.1 <- length(procs.list) < 2
  logic.2 <- ! names(procs.list) %in% c('direct', 'indirect')
  if ( logic.1 || logic.2 ) stop(err)

  direct.procs <- procs.list$direct
  indirect.procs <- procs.list$indirect

  # Checking for required procs entries
  direct.required <- c('sw.p', 'si', 'sf', 'reverse', 'offset')
  direct.procs <- .validate_param(direct.procs, direct.required)

  # Checking for required proc2s entries
  indirect.required <- c('sw.p', 'si', 'sf', 'reverse', 'offset')
  indirect.procs <- validate_param(indirect.procs, indirect.required)

  # Extracting parameters
  sw.p <- c(direct.procs$sw.p, indirect.procs$sw.p)
  si <- c(direct.procs$si, indirect.procs$si)
  sf <- c(direct.procs$sf, indirect.procs$sf)
  rv <- c(direct.procs$reverse, indirect.procs$reverse)
  ofs <- c(direct.procs$offset, indirect.procs$offset)

  n <- si[1]*si[2]

  # Doing some basic validation
  path <- .validate_pdata(path, number)

  # Checking which files actually exist
  components <- c('2rr', '2ri', '2ir', '2ii')
  all.paths <- file.path(path, components)
  existing.paths <- file.exists(all.paths)

  # Error out if there is no data found
  err <- 'At least one of either 2rr, 2ri, 2ir, or 2ii must exist to load data.' 
  if ( sum(existing.paths) == 0 ) stop(err)

  # Reading all available data
  f_read <- function(path) {

    if ( file.exists(path) ) {
      intensity <- safe_read(path, 'bin', size = 4, what = 'integer', n = n)

      # Checking load
      err <- 'Error reading processed files, file size does not match data.'
      if ( length(intensity) < n ) stop(err)
    
      intensity
    }
    else {
      0
    }
  }

  intensity <- lapply(all.paths, f_read)
  names(intensity) <- components

  # If either ir or ri is present, use cmplx2
  if ( any( c('2ri', '2ir') %in% components[existing.paths] ) ) {
    intensity <- cmplx2(rr = intensity[['2rr']], ri = intensity[['2ri']],
                        ir = intensity[['2ir']], ii = intensity[['2ii']])
  }
  # Otherwise, store 2rr and 2ii as simple r, i components
  else {
    intensity <- cmplx1(r = intensity[['2rr']], i = intensity[['2ii']])
  }

  # Formatting direct.shift 
  direct.shift <- seq(ofs[1], ofs[1] - sw.p[1]/sf[1], length.out = si[1])
  if (rv[1] == 'yes') direct.shift <- rev(direct.shift)

  # Formatting indirect.shift 
  indirect.shift <- seq(ofs[2] - sw.p[2]/sf[2], ofs[2], length.out = si[2])
  if (rv[2] == 'yes') indirect.shift <- rev(indirect.shift)

  # Combining output
  tibble(direct.shift = rep(direct.shift, si[2]), 
         indirect.shift = rep(indirect.shift, each = si[1]),
         intensity = intensity)

}




#==============================================================================>
#  Formatting and printing
#==============================================================================>



#' @export
format.NMRData2D <- function(x, ...) {
  d <- processed(x)
  components <- paste(colnames(as_tibble(d$intensity)), collapse = ', ') 
  direct.range <- range(d$direct.shift)
  direct.range <- sprintf("%.2f ppm to %.2f ppm direct", 
                          direct.range[1], direct.range[2])
  indirect.range <- range(d$indirect.shift)
  indirect.range <- sprintf("%.2f ppm to %.2f ppm indirect", 
                            indirect.range[1], indirect.range[2])
  sprintf('NMRData1D object (%s), %s, %s\n', 
          components, direct.range, indirect.range)
}

#' @export
print.NMRData2D <- function(x, ...) cat(format(x))

#' @export
setMethod("show", "NMRData2D", 
  function(object) cat(format(object))
  )

#' @export
summary.NMRData2D <- function(object, ...) {
  d <- object@processed
  summary(bind_cols(direct.shift = d$direct.shift,
                    indirect.shift = d$indirect.shift,
                    as_tibble(d$intensity)))
}



#' @export
is_vector_s3.NMRData2D <- function(x) FALSE

#' @export
type_sum.NMRData2D <- function(x) "NMRData2D"

#' @export
obj_sum.NMRData2D <- function(x) format(x)



#==============================================================================>
#  Plotting
#==============================================================================>

#------------------------------------------------------------------------------
#' Plot NMRData2D object
#' 
#' Convenience function that generates a plot of the spectral data.
#' 
#' @param x An NMRData2D object.
#' @param components Some combination of 'rr', 'ii', 'ir', and 'ri' in a single
#'                   string such as 'rr/ii' to include various real and
#'                   imaginary components. If more than one component is
#'                   selected, they are displayed in separate subplots.
#' 
#' @return A plot_ly plot.
#' 
#' @export
plot.NMRData2D <- function(x, components = 'rr') {

  legend.opts <- list(orientation = 'h', xanchor = "center", x = 0.5)

  f <- list(
    family = "Courier New, monospace",
    size = 16,
    color = "#7f7f7f"
  )

  xaxis <- list(
    title = "Indirect chemical shift (ppm)",
    titlefont = f
  )
  yaxis <- list(
    title = "Direct chemical shift (ppm)",
    titlefont = f,
    autorange = "reversed"
  )
  zaxis <- list(
    title = "Intensity",
    titlefont = f
  )

  scene <- list(
    xaxis = xaxis,
    yaxis = yaxis,
    zaxis = zaxis,
    camera=list(
      eye = list(x=0.5, y=0.5, z=0.5)
    )
  )

  #---------------------------------------

  d <- x@processed
  direct.shift <- d$direct.shift
  indirect.shift <- d$indirect.shift
  y.data <- d$intensity

  # Defining generic plot function
  f_init <- function(x, y, z, color, name, scene) {

    # Drawing separate lines for each indirect dimension
    groups <- unique(y)
    n <- length(groups)

    index <- y == groups[1]
    p <- plot_ly(x = x[index], y = y[index], z = z[index],
                 name = I(name), color = I(color), scene = scene,
                 type = 'scatter3d', mode = 'lines', legendgroup = name)


    # Looping over the rest
    for ( value in groups[-1] ) {
      index <- y == value
      p <- p %>%
        add_trace(x = x[index], y = y[index], z = z[index],
                  name = I(name), color = I(color), scene = scene,
                  type = 'scatter3d', mode = 'lines', legendgroup = name,
                  showlegend = FALSE)
    }

    p
  }

  # Initializing the plot list
  plots <- list()

  # Checking which components to plot
  rr <- grepl('rr', components)
  ri <- grepl('ri', components)
  ir <- grepl('ir', components)
  ii <- grepl('ii', components)

  # Plotting 
  if ( rr ) plots$rr <- f_init(direct.shift, indirect.shift, y.data$rr, 
                               'black', 'Real', 'scene1')
  if ( ri ) plots$ri <- f_init(direct.shift, indirect.shift, y.data$ri, 
                               'grey', 'Imaginary/Real', 'scene2')
  if ( ir ) plots$ir <- f_init(direct.shift, indirect.shift, y.data$ir, 
                               'grey', 'Real/Imaginary', 'scene3')
  if ( ii ) plots$ii <- f_init(direct.shift, indirect.shift, y.data$ii, 
                               'grey', 'Imaginary', 'scene4')

  if ( length(plots) == 0 )  return(NULL)
  else p <- subplot(plots, shareX = TRUE, shareY = TRUE, 
                    nrows = min(length(plots), 2))

# Tacking on scene settings
  args <- list(p)
  scenes <- c('scene1', 'scene2', 'scene3', 'scene4')
  for ( i in 1:length(plots) ) args[[scenes[i]]] <- scene

  do.call(layout, args) 
}

setMethod("plot", "NMRData2D", plot.NMRData2D)
