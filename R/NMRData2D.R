# Definition of a class structure for 2D NMR data.



#------------------------------------------------------------------------------
#' A class combining NMR data and scan parameters
#' 
#' An extension of the generic NMRData class to provide 2D-specific methods
#' such as constructors and plotting.
#' 
#' @rdname NMRData2D
#' @export
NMRData2D <- setClass("NMRData2D", contains = "NMRData")

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
#' Constructors for generating an NMRData2D object
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
#' @return An NMRData2D object.
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



#------------------------------------------------------------------------------
#' Convert NMRScaffold2D to NMRData2D
#' 
#' Generate lineshape from specified peaks.
#' 
#' @param object An NMRScaffold2D object
#' @param direct.shift Vector of chemical shift data in ppm.
#' @param indirect.shift Similar to direct.shift but for the indirect dimension.
#' @param direct.sf Sweep frequency (in MHz) -- needed to convert peak widths
#'                  from Hz to ppm. In most cases, it is recommended to set a
#'                  single default value using nmroptions$direct$sf = ..., but
#'                  an override can be provided here.
#' @param indirect.sf Similar to direct.sf but for the indirect dimension.
#' 
#' @return An NMRData2D object.
#' 
#' @export
nmrdata_2d_from_scaffold <- function(
  object, direct.shift = NULL, indirect.shift = NULL, 
  direct.sf = nmroptions$direct$sf, indirect.sf = nmroptions$indirect$sf) {

  if ( is.null(direct.shift) ) {
    positions <- peaks(direct(object))$position
    direct.shift <- seq(min(positions) - 0.2, max(positions) + 0.2,
                        length.out = 100)
  }

  if ( is.null(indirect.shift) ) {
    positions <- peaks(indirect(object))$position
    indirect.shift <- seq(min(positions) - 0.2, max(positions) + 0.2,
                        length.out = 100)
  }

  # Generate grid
  logic.1 <- length(direct.shift) == length(unique(direct.shift))
  logic.2 <- length(indirect.shift) == length(unique(indirect.shift))
  if ( logic.1 && logic.2 ) {
    d <- expand.grid(x1 = direct.shift, x2 = indirect.shift)
    direct.shift <- d$x1
    indirect.shift <- d$x2
  }

  # Using the procs file to load the processed data
  processed <- values(object, direct.shift = direct.shift,
                      indirect.shift = indirect.shift,
                      direct.sf = direct.sf, indirect.sf = indirect.sf,
                      use.cmplx1 = TRUE)

  # acqus just contains the sf
  acqus <- list(direct = list(sf = direct.sf), indirect = list(sf = indirect.sf))

  # Returning class object
  new("NMRData2D", processed = processed, acqus = acqus)
}



#==============================================================================>
#  Processing
#==============================================================================>



#------------------------------------------------------------------------
# Extracting lineshape values
#' @rdname values
#' @export
setMethod("values", "NMRData2D", 
  function(object, domain = "rr/ri/ir/ii", use.cmplx1 = FALSE) {
    
    d <- object@processed
    d$intensity <- domain(d$intensity, domain, use.cmplx1)

    d
  })





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
  sprintf('NMRData2D object (%s), %s, %s\n', 
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
range.NMRData2D <- function(object, ...) {
  d <- object@processed
  list(direct.shift = range(d$direct.shift),
       indirect.shift = range(d$indirect.shift),
       intensity = range(Re(d$intensity)))
}

#' @export
is_vector_s3.NMRData2D <- function(x) FALSE

#' @export
type_sum.NMRData2D <- function(x) "NMRData2D"

#' @export
obj_sum.NMRData2D <- function(x) format(x)




#==============================================================================>
# Access to direct and indirect components
#==============================================================================>



#------------------------------------------------------------------------------
# Projection

#---------------------------------------
#' Get 1D projection of dimension
#' 
#' Warning: this is currently a convenience method to facilitate internal data
#' manipulation, so the conversion may drop important data. When applied to a
#' an NMRScaffold2D object, non-selected components are simply dropped, and
#' when applied to an NMRData2D object, the non-selected are summed together to
#' give a 1D projection. No normalization is performed so non-uniform sampling
#' may result in distorted results. The functions direct(object) and
#' indirect(object) are provided as shortcuts to projection(object, "direct")
#' and projection(object, "indirect").
#' 
#' @param object An NMRData2D or NMRScaffold2D object.
#' @param dimension Currently limited to either "direct" or "indirect"
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name project
setGeneric("projection", 
  function(object, dimension, ...) standardGeneric("projection")
)

#' @rdname projection
#' @export
setMethod("projection", "NMRData2D",
  function(object, dimension) {

  intensity <- as_tibble(object@processed$intensity)

  # Double checking if dimension exists
  if ( dimension == "direct" ) {
    r.col <- pull(intensity, "rr")
    i.col <- pull(intensity, "ir")
  } else if ( dimension == "indirect" ) {
    r.col <- pull(intensity, "rr")
    i.col <- pull(intensity, "ri")
  } else {
    err <- sprintf('Dimension "%s" not found', dimension)
    stop(err)
  }

  x <- paste(dimension, "shift", sep = ".")

  # Parameters (revisit and cross-reference NMRData)
  procs <- list(direct = list())#, indirect = list())
  procs$direct <- object@procs[[dimension]]

  acqus <- list(direct = list())#, indirect = list())
  acqus$direct <- object@acqus[[dimension]]

  d <- object@processed
  d$r <- r.col
  d$i <- i.col

  # Handling processed data
  d <- d %>%
    select(-intensity) %>%
    group_by(across(x)) %>%
    summarize(r = sum(r), i = sum(i)) %>%
    ungroup()

  d <- tibble(direct.shift = d[[x]], intensity = cmplx1(r = d$r, i = d$i))

  new("NMRData1D", processed = d, procs = procs, acqus = acqus)
})



#' @rdname projection
#' @export
setGeneric("direct", 
  function(object) standardGeneric("direct")
)

#' @rdname projection
#' @export
setMethod("direct", "NMRData2D",
  function(object) {
  
  projection(object, "direct")
})



#' @rdname projection
#' @export
setGeneric("indirect", 
  function(object) standardGeneric("indirect")
)

#' @rdname projection
#' @export
setMethod("indirect", "NMRData2D",
  function(object) {
  
  projection(object, "indirect")
})



#==============================================================================>
#  Plotting
#==============================================================================>



#------------------------------------------------------------------------------
#' Plot NMRData2D object
#' 
#' Convenience function that generates a plot of the spectral data.
#' 
#' @param x An NMRData2D object.
#' @param components One of 'rr', 'ii', 'ir', or 'ri' to specify various real
#'                   and imaginary components. 3D plot subplots are not
#'                   currently supported.
#' 
#' @return A plot_ly plot.
#' 
#' @export
plot.NMRData2D <- function(x, components = 'rr') {

  components <- tolower(components)
  valid.components <- c('rr', 'ri', 'ir', 'ii')
  err <- '"components" must be one of "rr", "ri", "ir", or "ii".'
  if (! components %in% valid.components ) stop(err)

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
    title = "Indirect chemical shift (ppm)",
    titlefont = f,
    autorange = "reversed"
  )
  zaxis <- list(
    title = "Intensity",
    titlefont = f
  )

  scene <- list(
    legend = legend.opts,
    xaxis = xaxis,
    yaxis = yaxis,
    zaxis = zaxis,
    camera=list(
      eye = list(x=0, y=2, z=1.25)
    )
  )

  #---------------------------------------

  d <- x@processed
  direct.shift <- d$direct.shift
  indirect.shift <- d$indirect.shift
  y.data <- d$intensity

  # Defining generic plot function
  f_init <- function(x, y, z, color, name) {

    # Drawing separate lines for each indirect dimension
    groups <- unique(y)
    n <- length(groups)

    index <- y == groups[1]
    p <- plot_ly(x = x[index], y = y[index], z = z[index],
                 name = I(name), color = I(color),
                 type = 'scatter3d', mode = 'lines', legendgroup = name)


    # Looping over the rest
    for ( value in groups[-1] ) {
      index <- y == value
      p <- p %>%
        add_trace(x = x[index], y = y[index], z = z[index],
                  name = I(name), color = I(color),
                  type = 'scatter3d', mode = 'lines', legendgroup = name,
                  showlegend = FALSE)
    }

    p %>%
      layout(scene = scene)
  }

  # Initializing the plot list
  plots <- list()

  # Checking which components to plot
  rr <- grepl('rr', components)
  ri <- grepl('ri', components)
  ir <- grepl('ir', components)
  ii <- grepl('ii', components)

  # Plotting 
  x <- direct.shift
  y <- indirect.shift
  p <- NULL
  if ( rr ) p <- f_init(x, y, y.data$rr, 'black', 'Real')
  if ( ri ) p <- f_init(x, y, y.data$ri, 'black', 'Real/Imaginary')
  if ( ir ) p <- f_init(x, y, y.data$ir, 'black', 'Imaginary/Real')
  if ( ii ) p <- f_init(x, y, y.data$ii, 'black', 'Imaginary')

  p
}

setMethod("plot", "NMRData2D", plot.NMRData2D)
