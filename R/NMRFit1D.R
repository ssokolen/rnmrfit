# Definition of a class structure for lineshape fitting.



#==============================================================================>
#  NMRFit1D
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR fit.
#' 
#' This class is used to extend an NMRMixture1D object with NMRData1D while
#' also defining baseline and phase correction terms. There is just one primary
#' method associated with this class: fit(). fit() takes input data and peak
#' definitions, applies a nonlinear least squares fit, and generates best-fit
#' peak parameters, overwriting the initial values. Since the fit process is
#' destructive, the nmrfit_1d() function used to initialize a fit object has an
#' option to delay the fit, allowing pre-fit and post-fit objects to be saved
#' as different variables.
#' 
#' @slot nmrdata An NMRData1D object used to fit the peaks.
#' @slot knots A vector of chemical shifts corresponding to the internal knots
#'             of the baseline term (two boundary knots are always included).
#' @slot baseline A vectors of complex numeric values representing the baseline.
#'                The actual baseline is generated from a basis spline so these
#'                values are control points that are roughly proportional to the
#'                intensity values of the baseline based on the chemical shift
#'                locations of the knots. The default baseline order and the
#'                number of internal knots can be set using
#'                nmrpotions$baseline$order = 3, nmroptions$baseline$n.knots =
#'                0, with the baseline vector length equal to the order +
#'                n.knots + 1.
#' @slot phase A vector of phase correction terms, corresponding to either 0 or
#'             1st order phase correction. The default phase correction order
#'             can be set using nmroptions$phase$order = 0. Note that the 0
#'             order term is always calculated with respect to 0 ppm and the
#'             first order term has units of degrees per ppm.
#' @slot bounds A lower and upper bounds on baseline and phase terms. Both lower
#'              and upper bounds are lists containing "baseline", and "phase"
#'              elements. A "peaks" element is also generated dynamically from
#'              the species list when using the bounds() getter function. Unlike
#'              peaks, baseline and phase only have a single lower and upper
#'              bound, representing the overall minimum/maximum baseline/phase
#'              value at any point in the spectrum.
#' 
#' @name NMRFit1D-class
#' @export
NMRFit1D <- setClass("NMRFit1D",
  contains = 'NMRMixture1D',
  slots = c(
    nmrdata = 'NMRData1D',
    knots = 'numeric',
    baseline = 'complex',
    phase = 'numeric',
    lower.bounds = 'list',
    upper.bounds = 'list',
    time = 'numeric',
    status = 'character',
    code = 'numeric'
  ),
  prototype = prototype(
    knots = numeric(0), 
    baseline = complex(re = rep(0, 4), im = rep(0, 4)),
    phase = c(0),
    lower.bounds = NULL,
    upper.bounds = NULL,
    time = c(0),
    status = 'no fit attempted',
    code = c(0)
  )
)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRFit1D object
#' 
#' Generates an NMRFit1D object based on an NMRMixture1D or a list of
#' NMRSpecies1D objects. See ?nmrmixture_1d and ?nmrspecies_1d for more details
#' about this conversion. It is important to note that all NMR objects defined
#' in the species list must fall inside the provided data (although partial
#' fitting may be implemented in the future).
#' 
#' @param species A list of NMRSpecies1D objects or other objects that can be
#'                converted to NMRSpecies1D objects. See ?nmrresonance_1d and
#'                ?nmrspecies_1d for more details about this conversion. If list
#'                elements are named, these names will be use to replace species
#'                ids.
#' @param nmrdata An NMRData1D object used to fit the supplied peaks.
#'                automatically generated from the resonance names.
#' @param baseline.order An integer specifying the order of the baseline spline
#'                       function. Note the the internal B-spline implementation
#'                       requires an order of 1 or greater. You can use an order
#'                       of -1 to disable baseline correction.
#' @param n.knots An integer specifying the number of internal knots to use. The
#'                specific position of these knots can be modified later using
#'                knots() function. To modify initial values of the baseline,
#'                use baseline() function.
#' @param phase.order An integer specifying the order of the phase correction
#'                    polynomial -- may be 0 or 1 (-1 to disable). in practice,
#'                    but a higher order correction is possible.
#' @param delay.fit FALSE to immediately run least squares optimization after
#'                  the NMRFit1D object is initialized, TRUE to skip the
#'                  optimization, enabling more customization. The fit can be
#'                  run manually using fit().
#' @param init An initialization, function that takes an NMRFit1D object and
#'             returns a modified NMRFit1D object. Use the "identity" function
#'             to override the default initialization in the
#'             nmroptions$fit$init. Note that all arguments provided to the
#'             fit() function are also passed on the init() function.
#' @param opts A list of NLOPT fit options to override the default options in
#'             the nmroptions$fit$opts.
#' @param ... Options passed to nmrmixture_1d if conversion has to be performed.
#'            See ?nmrmixture_1d for more details.
#' 
#' @return An NMRFit1D object.
#' 
#' @export
nmrfit_1d <- function(
  species, nmrdata, baseline.order = nmroptions$baseline$order,
  n.knots = nmroptions$baseline$n.knots, 
  phase.order = nmroptions$phase$order, 
  delay.fit = FALSE, init = nmroptions$fit$init, 
  opts = nmroptions$fit$opts, ...) {

  #---------------------------------------
  # Generating list of species 

  # Generate an NMRMixture object that will be used as a template for the fit
  mixture <- nmrmixture_1d(species, ...)

  #---------------------------------------
  # Checking nmrdata

  check_conformity(mixture, nmrdata)

  #---------------------------------------
  # Baseline and phase

  # The initial value for the baseline is just the median of nmrdata intensity
  if ( baseline.order == -1 ) {
    baseline <- complex(0)
    knots <- numeric(0)
  }
  else {
    n <- baseline.order + n.knots + 1
    baseline <- complex(re = rep(median(Re(nmrdata@processed$intensity)), n),
                        im = rep(median(Im(nmrdata@processed$intensity)), n))

    # Knots are initialized to fall evenly between the chemical shift data
    if ( n.knots > 0 ) {
      knots <- seq(0, 1, length.out = n.knots + 2)
      knots <- knots[2:(n.knots + 1)]
    } else {
      knots <- numeric(0)
    }
  }

  # The initial value for the phase is always 0
  if ( phase.order == -1 ) phase <- numeric(0)
  else if ( phase.order > 1 ) {
    err <- "Phase order must be 0 or 1 (or -1 to disable phase correction)"
    stop(err)
  }
  else phase <- rep(0, phase.order + 1)

  # Initializing bounds
  lower_bounds = list(baseline = complex(0), phase = numeric(0))
  upper_bounds = list(baseline = complex(0), phase = numeric(0))

  #---------------------------------------
  # Resulting fit object
  out <- new('NMRFit1D', mixture, nmrdata = nmrdata,
                         knots = knots, baseline = baseline, phase = phase,
                         lower.bounds = lower_bounds, 
                         upper.bounds = upper_bounds)

  # If the fit is delayed, then return current object, otherwise run fit first
  if ( delay.fit ) out
  else fit(out, init = init, opts = opts)
}



#==============================================================================>
#  The main fit function (wrapping Rcpp code)
#==============================================================================>



#------------------------------------------------------------------------------
#' Fit NMR data to peaks.
#' 
#' Given an NMRFit1D object, this function performs a least squares fit on the
#' data and updates the peak parameters. 
#' 
#' @param object An NMRFit1D object.
#' @param init An initialization function that takes NMRFit1D object and returns
#'             a modified NMRFit1D object. Use the "identity" function to
#'             override the default initialization in the
#'             nmroptions$fit$init. Note that all arguments provided to the
#'             fit() function are also passed on the init() function.
#' @param opts A list of NLOPT fit options to override the default options in
#'             the nmroptionsd$fit$opts.
#' 
#' @name fit
#' @export
#' @useDynLib rnmrfit fit_1d_wrapper
setGeneric("fit", 
  function(object, ...) {
    standardGeneric("fit")
})

#' @rdname fit
#' @export
setMethod("fit", "NMRFit1D",
  function(object, init = nmroptions$fit$init, opts = nmroptions$fit$opts) {

  # First, run the initialization
  object <- init(object)

  # Ensuring consistent order
  d <- processed(object@nmrdata)
  d <- d[order(d$direct.shift), ]

  # Unpacking some of the NMR data
  x <- d$direct.shift
  x.range <- range(x)
  x.span <- x.range[2] - x.range[1]
  y <- d$intensity
  y.range <- range(Re(y))

  sf <- object@sf

  # Normalizing data
  x <- (x - x.range[1])/x.span
  y <- complex(re = Re(y), im = Im(y))
  y <- y/y.range[2]

  # Scaling and unpacking all parameters
  peaks <- peaks(object)
  n.peaks <- nrow(peaks)
  data.columns <- c('position', 'width', 'height', 'fraction.gauss')

  peaks <- list(par = peaks[, data.columns], 
                lb = bounds(object)$lower$peaks[, data.columns], 
                ub = bounds(object)$upper$peaks[, data.columns])

  for (name in names(peaks)) {
    peaks[[name]]$position <- (peaks[[name]]$position - x.range[1])/x.span
    peaks[[name]]$width <- peaks[[name]]$width/sf/x.span
    peaks[[name]]$height <- peaks[[name]]$height/y.range[2]

    peaks[[name]] <- as.vector(t(as.matrix(peaks[[name]])))
  }

  #---------------------------------------
  # Scaling and unpacking baseline/phase terms

  # Tacking out bounds to internal knots
  knots <- sort(c(knots(object), c(0, 1)))
  
  n.baseline <- length(baseline(object))

  # y-scaling performed below for convenience
  baseline <- list(par = baseline(object), 
                   lb = rep(bounds(object)$lower$baseline, n.baseline),
                   ub = rep(bounds(object)$upper$baseline, n.baseline))

  # 1st order phase coefficients must be adapted to the local scale
  # (although a 0 order correction remains constant)
  phase <- phase(object)
  n.phase <- length(phase)
  if ( n.phase == 2 ) {
    phase[1] <- phase[1] + phase[2]*x.range[1]
    phase[2] <- (phase[2] - phase[1])*x.span
  }

  # Adding simple bounds that are expanded to constraints for 1st order
  # phase correction
  phase <- list(par = phase, 
                lb = rep(bounds(object)$lower$phase, n.phase),
                ub = rep(bounds(object)$upper$phase, n.phase))

  par <- list(par = NA, lb = NA, ub = NA)
  for (name in names(par)) {

    # At the end of the day, all parameters must be flattened into a vector
    par[[name]] <- c(
      peaks[[name]], 
      Re(baseline[[name]])/y.range[2], 
      Im(baseline[[name]])/y.range[2], 
      phase[[name]]
    )
  }

  #---------------------------------------
  # Checking opts
  algorithm <- opts$algorithm
  err <- '"algorithm" option must be one of "SLSQP" or "ISRES"'
  if (! algorithm %in% c("SLSQP", "ISRES") ) stop(err)

  if ( algorithm == "SLSQP" ) algorithm <- 0
  else algorithm <- 1

  padding <- opts$padding/x.span
  err <- '"padding" option must 0 or greater'
  if ( padding < 0 ) stop(err)

  xtol_rel <- opts$xtol_rel
  err <- '"xtol_rel" option must be between 0 and 1'
  if ( (xtol_rel <= 0) | (xtol_rel >= 1) ) stop(err)

  maxtime <- opts$maxtime
  err <- '"maxtime" option must be 0 (disabled) or greater'
  if ( maxtime < 0 ) stop(err)

  #---------------------------------------
  # Generating constraint lists

  constraints <- parse_constraints(object, x.span)

  eq.constraints <- constraints[[1]]
  ineq.constraints <- constraints[[2]]

  # Adding padding constraints
  if ( padding > 0 ) {
    peaks <- peaks(object)
    n <- nrow(peaks)

    pos.left <- peaks$position[1:(n-1)]
    pos.right <- peaks$position[2:n]

    logic <- pos.left < pos.right
    pos.left <- ifelse(logic, 1:(n-1), -(1:(n-1)))
    pos.right <- ifelse(!logic, 2:n, -(2:n))

    pad.constraints <- map2(pos.left, pos.right, ~ c(0, -padding, .x, .y))
    ineq.constraints <- c(ineq.constraints, pad.constraints)
  }

  #---------------------------------------
  # Performing the fit

  # Flattening constraints
  eq.constraints <- unlist(lapply(eq.constraints, function (x) c(x, NaN)))
  eq.constraints <- eq.constraints[-length(eq.constraints)]

  ineq.constraints <- unlist(lapply(ineq.constraints, function (x) c(x, NaN)))
  ineq.constraints <- ineq.constraints[-length(ineq.constraints)]

  start.time <- proc.time()
  out <- .Call("fit_1d_wrapper", 
    x = as.double(x), 
    y = as.double(c(Re(y), Im(y))), 
    knots = as.double(knots),
    p = as.double(par$par), 
    lb = as.double(par$lb), 
    ub = as.double(par$ub), 
    n = as.integer(length(x)),
    nl = as.integer(n.peaks*4),
    nb = as.integer(n.baseline),
    np = as.integer(n.phase),
    nk = as.integer(length(knots)),
    eq = as.double(eq.constraints),
    iq = as.double(ineq.constraints),
    neq = as.integer(length(eq.constraints)),
    niq = as.integer(length(ineq.constraints)),
    alg = as.integer(algorithm),
    xtr = as.double(xtol_rel),
    mxt = as.double(maxtime),
    out = as.integer(0)
  )
  object@time <- as.numeric(proc.time() - start.time)[3]

  # Defining error codes
  positive_codes <- c(
    "success", "stop value reached", "ftol value reached", "xtol reached",
    "max evaluations reached", "max time reached"
  )

  negative_codes <- c(
    "unspecified failure", "invalid arguments", "out of memory",
    "significant roundoff errors", "forced stop"
  )

  code <- out[[2]]
  print(code)
  if ( (code > 0) & (code <= length(positive_codes)) ) {
    object@status <- positive_codes[code]
    object@code <- 0
  } else if ( (code < 0) & (-code <= length(negative_codes)) ) {
    object@status <- negative_codes[-code]
    object@code <- -code
  } else {
    object@status <- "error communicating with rust"
    object@code <- 99
  }

  # Issuing warning messages as required
  if ( code == 0 ) {
    wrn <- "Error communicating with rust -- solution may be suboptimal."
    warning(wrn)
  } else if ( code > 4 ) {
    wrn <- "Fit terminated prematurely -- %s -- solution may be suboptimal"
    wrn <- sprintf(wrn, object@status)
    warning(wrn)
  } else if ( code < 0 ) {
    wrn <- "Fit failed -- %s"
    wrn <- sprintf(wrn, object@status)
    warning(wrn)
  }

  par$par <- out[[1]]

  #---------------------------------------
  # Unpacking and rescaling parameters
  
  # Starting with peaks
  peaks <- peaks(object)
  new.peaks <- matrix(par$par[1:(n.peaks*4)], ncol = 4, byrow = TRUE)
  peaks[, data.columns] <- new.peaks

  peaks$position <- peaks$position*x.span + x.range[1]
  peaks$width <- peaks$width*sf*x.span
  peaks$height <- peaks$height*y.range[2]

  peaks(object) <- peaks

  # Then baseline
  index.re <- (n.peaks*4 + 1):(n.peaks*4 + n.baseline)
  index.im <- index.re + n.baseline
  if ( n.baseline > 0 ) {
    object@baseline <- complex(re = par$par[index.re]*y.range[2],
                               im = par$par[index.im]*y.range[2])
  }

  # And phase
  index <- (n.peaks*4 + n.baseline*2 + 1):(n.peaks*4 + n.baseline*2 + n.phase)
  if ( n.phase > 0 ) {
    phase <- par$par[index]

    # 1st order phase coefficients must be adapted back to the global scale
    if ( n.phase == 2 ) {
      phase.left <- phase[1]
      phase.right <- phase[1] + phase[2]
      phase[2] <- (phase.right - phase.left)/x.span
      phase[1] <- phase.left - phase[2]*x.range[1]
    }

    object@phase <- phase
  }

  object
})




#==============================================================================>
#  Display function
#==============================================================================>



#------------------------------------------------------------------------------
#' Display NMRFit1D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRFit1D", 
  function(object) {

    # Generating compiled data frames
    peaks <- peaks(object)
    baseline <- baseline(object)
    knots <- knots(object)
    phase <- phase(object)
    lower <- lower_bounds(object)
    upper <- upper_bounds(object)
    couplings <- couplings(object)

    cat('An object of NMRFit1D class\n\n')

    # Peaks
    cat('Peaks:\n\n')
    print(peaks)
    cat('\n')

    # Baseline
    cat('Baseline correction:\n\n')

    cat('Real baseline: ')
    if ( length(baseline) > 0) {cat('\n'); print(round(Re(baseline), 4))}
    else cat('None\n')

    cat('Internal knots: ')
    if ( length(knots) > 0) {cat('\n'); print(round(knots, 4))}
    else cat('None\n')
    cat('\n')

    # Phase
    cat('Phase correction:\n\n')

    if ( length(phase) > 0) print(round(phase, 4))
    else cat('None\n')
    cat('\n')

    # Bounds
    columns <- c('position', 'width', 'height', 'fraction.gauss')
    for ( column in columns ) {
      range <- paste('(', lower[, column], ', ', 
                          upper[, column], ')', sep = '')
      peaks[, column] <- range
    }

    cat('Bounds (lower, upper):\n\n')
    print(peaks)
    cat('\n')   

    # Baseline bounds
    cat('Baseline and phase correction bounds:\n\n')

    cat('Real baseline: ')
    if ( length(object@lower.bounds$baseline) > 0) {
      lower <- round(Re(object@lower.bounds$baseline), 4)
      upper <- round(Re(object@upper.bounds$baseline), 4)
      
      range <- paste('(', lower, ', ', upper, ')\n', sep = '')
      cat(range)
    }
    else cat('None\n')

    # Phase bounds
    cat('Phase (radians): ')

    if ( length(object@lower.bounds$phase) > 0)  {
      lower <- round(object@lower.bounds$phase, 4)
      upper <- round(object@upper.bounds$phase, 4)
      
      range <- paste('(', lower, ', ', upper, ')\n', sep = '')
      cat(range)
    }
    else cat('None\n')
    cat('\n')

    # Couplings
    if ( nrow(couplings) > 0 ) {
      cat('Couplings:\n\n')
      print(couplings)
      cat('\n')
    }
    else {
      cat('No couplings defined.\n')
    }

  })



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
# Bounds

#' @rdname bounds
#' @export
setMethod("bounds", "NMRFit1D", 
  function(object) {
    # Extracting peak bounds from species using the NMRMixture1D signature
    bounds <- selectMethod("bounds", signature="NMRScaffold1D")(object)

    # Baseline and phase
    lower.baseline <- object@lower.bounds$baseline
    lower.baseline <- ifelse(length(lower.baseline) == 0, 
                             complex(re = -Inf, im = -Inf), lower.baseline)
    upper.baseline <- object@upper.bounds$baseline
    upper.baseline <- ifelse(length(upper.baseline) == 0, 
                             complex(re = Inf, im = Inf), upper.baseline)   

    lower.phase <- object@lower.bounds$phase
    lower.phase <- ifelse(length(lower.phase) == 0, -Inf, lower.phase)   
    upper.phase <- object@upper.bounds$phase
    upper.phase <- ifelse(length(upper.phase) == 0, Inf, upper.phase)   

    # Outputting
    list(lower = list(peaks = bounds$lower, baseline = lower.baseline, 
                      phase = lower.phase), 
         upper = list(peaks = bounds$upper, baseline = upper.baseline,
                      phase = upper.phase))
  })



#------------------------------------------------------------------------------
# Baseline

#---------------------------------------
#' Get object baseline 
#' 
#' Generic convenience method to access the baseline definition of an
#' NMRFit1D object.
#' 
#' @param object An NMRFit1D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name baseline
#' @export
setGeneric("baseline", 
  function(object, ...) standardGeneric("baseline")
)

#' @rdname baseline
#' @export
setMethod("baseline", "NMRFit1D", 
  function(object) object@baseline
)



#------------------------------------------------------------------------------
# Knots

#---------------------------------------
#' Get object baseline knots 
#' 
#' Generic convenience method to access the baseline knots definition of an
#' NMRFit1D object.
#' 
#' @param object An NMRFit1D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name knots
#' @export
setGeneric("knots", 
  function(object, ...) standardGeneric("knots")
)

#' @rdname knots
#' @export
setMethod("knots", "NMRFit1D", 
  function(object) object@knots
)



#------------------------------------------------------------------------------
# Phase

#---------------------------------------
#' Get object phase 
#' 
#' Generic convenience method to access the phase definition of an
#' NMRFit1D object.
#' 
#' @param object An NMRFit1D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name phase
#' @export
setGeneric("phase", 
  function(object, ...) standardGeneric("phase")
)

#' @rdname phase
#' @export
setMethod("phase", "NMRFit1D", 
  function(object) object@phase
)
