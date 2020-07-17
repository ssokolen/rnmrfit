# Definition of a class structure for lineshape fitting.



#==============================================================================>
#  NMRFit2D
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR fit.
#' 
#' This class is used to extend an NMRMixture2D object with NMRData2D while
#' also defining baseline and phase correction terms. There is just one primary
#' method associated with this class: fit(). fit() takes input data and peak
#' definitions, applies a nonlinear least squares fit, and generates best-fit
#' peak parameters, overwriting the initial values. Since the fit process is
#' destructive, the nmrfit_2d() function used to initialize a fit object has an
#' option to delay the fit, allowing pre-fit and post-fit objects to be saved
#' as different variables.
#' 
#' @slot nmrdata An NMRData2D object used to fit the peaks.
#' @slot knots A vector of chemical shifts corresponding to the internal knots
#'             of the baseline term (two boundary knots are always included).
#' @slot baseline A vectors of complex numeric values representing the baseline.
#'                The actual baseline is generated from a basis spline so these
#'                values are control points that are roughly proportional to the
#'                intensity values of the baseline based on the chemical shift
#'                locations of the knots. The default baseline order and the
#'                number of internal knots can be set using
#'                nmrpotions$direct$baseline = list(order = 3, n.knots = 0),
#'                with the baseline vector length equal to the order + n.knots +
#'                1.
#' @slot phase A vector of phase correction terms, corresponding to either 0 or
#'             1st order phase correction. The default phase correction order
#'             can be set using nmroptions$direct$phase = list(order=1). Note
#'             that the 0 order term is always calculated with respect to 0 ppm
#'             and the first order term has units of degrees per ppm.
#' @slot bounds A lower and upper bounds on baseline and phase terms. Both lower
#'              and upper bounds are lists containing "baseline", and "phase"
#'              elements. A "peaks" element is also generated dynamically from
#'              the species list when using the bounds() getter function. Unlike
#'              peaks, baseline and phase only have a single lower and upper
#'              bound, representing the overall minimum/maximum baseline/phase
#'              value at any point in the spectrum.
#' 
#' @name NMRFit2D-class
#' @export
NMRFit2D <- setClass("NMRFit2D",
  contains = 'NMRMixture2D',
  slots = c(
    nmrdata = 'NMRData2D',
    knots = 'numeric',
    baseline = 'vctrs_cmplx2',
    phase = 'numeric',
    lower.bounds = 'list',
    upper.bounds = 'list',
    time = 'numeric'
  ),
  prototype = prototype(
    knots = numeric(0), 
    baseline = cmplx2(rr = rep(0, 4), ri = rep(0, 4),
                      ir = rep(0, 4), ii = rep(0, 4)),
    phase = c(0),
    lower.bounds = NULL, 
    upper.bounds = NULL,
    time = c(0)
  )
)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRFit2D object
#' 
#' Generates an NMRFit2D object based on a list of NMRSpecies2D objects or
#' other objects that can be converted to NMRSpecies2D objects. See
#' ?nmrresonance_2d and ?nmrspecies_2d for more details about this conversion.
#' Apart from providing peak definitions (via species argument) and data (via
#' nmrdata argument), the main fit decisions are related to baseline and phase
#' correction, as well as how to deal with peaks that are defined outside the
#' range of the data. This latter decision is broken down into two arguments:
#' exclusion.level and exclusion.notification. The exclusion.level parameter
#' determines which part of the overall species to exclude if any of its peaks
#' fall outside the data range: either 'species' for whole species, 'resonance'
#' for just a subset of the species and 'peak' to ignore resonance/species
#' blocks and exclude by specific peak alone. The exclusion.notification
#' parameter determines how to report when peaks are found to be outside the
#' data range: either 'none' to give no notice, 'message' to issue a message,
#' 'warning' to issue a warning, or 'stop' to issue an error.
#' 
#' @param species A list of NMRSpecies2D objects or other objects that can be
#'                converted to NMRSpecies2D objects. See ?nmrresonance_2d and
#'                ?nmrspecies_2d for more details about this conversion. If list
#'                elements are named, these names will be use to replace
#'                species ids.
#' @param nmrdata An NMRData2D object used to fit the supplied peaks.
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
#'                    polynomial -- may be 0 or 1 (-1 to disable).
#' @param delay.fit FALSE to immediately run least squares optimization after
#'                  the NMRFit2D object is initialized, TRUE to skip the
#'                  optimization, enabling more customization. The fit can be
#'                  run manually using fit().
#' @param direct.sf Sweep frequency (MHz) in the direct dimension -- needed to
#'                  convert coupling constants from Hz to ppm. In most cases, it
#'                  is recommended to set a single default value using
#'                  nmroptions$direct$sf  = ..., but an override can be
#'                  provided here.
#' @param indirect.sf Sweep frequency (MHz) in the indirect dimension -- needed
#'                    to convert coupling constants from Hz to ppm. In most
#'                    cases, it is recommended to set a single default value
#'                    using nmroptions$indirect$sf  = ..., but an override
#'                    can be provided here.
#' @param init An initialization, function that takes an NMRFit2D object and
#'             returns a modified NMRFit2D object. Use the "identity" function
#'             to override the default initialization in the
#'             nmroptions$fit$init. Note that all arguments provided to the
#'             fit() function are also passed on the init() function.
#' @param opts A list of NLOPT fit options to override the default options in
#'             the nmroptions$fit$opts.
#' @param ... Options passed to nmrspecies_2d if conversion has to be performed.
#'            See ?nmrspecies_2d for more details.
#' 
#' @return An NMRFit2D object.
#' 
#' @export
nmrfit_2d <- function(
  species, nmrdata, baseline.order = nmroptions$baseline$order,
  n.knots = nmroptions$baseline$n.knots, 
  phase.order = nmroptions$phase$order, 
  delay.fit = FALSE, direct.sf = nmroptions$direct$sf, 
  indirect.sf = nmroptions$indirect$sf, 
  init = nmroptions$fit$init, opts = nmroptions$fit$opts, ...) {

  # Checking to make sure that sweep frequency is defined
  err <- '"direct.sf" must be provided as input or set using nmroptions'
  if ( is.null(direct.sf) ) stop(err)

  # Checking to make sure that sweep frequency is defined
  err <- '"indirect.sf" must be provided as input or set using nmroptions'
  if ( is.null(indirect.sf) ) stop(err)

  #---------------------------------------
  # Generating list of species 

  # Generate an NMRMixture object that will be used as a template for the fit
  mixture <- nmrmixture_2d(species, ...)

  #---------------------------------------
  # Checking nmrdata

  check_conformity(mixture, nmrdata)

  #---------------------------------------
  # Baseline and phase

  # The initial value for the baseline is just the median of nmrdata intensity
  if ( baseline.order == -1 ) {
    baseline <- cmplx2(0)
    knots <- numeric(0)
  }
  else {
    n <- baseline.order + n.knots + 1
    baseline <- cmplx2(rr = rep(median(nmrdata@processed$intensity$rr), n),
                       ri = rep(median(nmrdata@processed$intensity$ri), n),
                       ir = rep(median(nmrdata@processed$intensity$ir), n),
                       ii = rep(median(nmrdata@processed$intensity$ii), n))

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
  else if ( phase.order == 0 ) phase <- rep(0, 2)
  else if ( phase.order == 1 ) phase <- rep(0, 3)
  else {
    err <- "Phase order must be 0 or 1 (or -1 to disable phase correction)"
    stop(err)
  }

  # Initializing bounds
  lower_bounds <- list(baseline = NULL, phase = numeric(0))
  upper_bounds <- list(baseline = NULL, phase = numeric(0))

  #---------------------------------------
  # Resulting fit object
  out <- new('NMRFit2D', mixture, nmrdata = nmrdata,
                         knots = knots, baseline = baseline, phase = phase,
                         lower.bounds = lower_bounds,
                         upper.bounds = upper_bounds)

  # If the fit is delayed, then return current object, otherwise run fit first
  if ( delay.fit ) out
  else fit(out, direct.sf = direct.sf, indirect.sf = indirect.sf, 
           init = init, opts = opts)
}



#==============================================================================>
#  The main fit function (wrapping Rcpp code)
#==============================================================================>



#' @rdname fit
#' @export
setMethod("fit", "NMRFit2D",
  function(object, direct.sf = nmroptions$direct$sf, 
           indirect.sf = nmroptions$indirect$sf, init = nmroptions$fit$init, 
           opts = nmroptions$fit$opts) {

  # First, run the initialization
  object <- init(object, sf = sf, init = init, opts = opts)

  # Ensuring consistent order
  d <- processed(object@nmrdata)
  d <- d[order(d$direct.shift, d$indirect.shift), ]

  # Unpacking some of the NMR data
  x1 <- d$direct.shift
  x2 <- d$indirect.shift
  x1.range <- range(x1)
  x2.range <- range(x2)
  x1.span <- x1.range[2] - x1.range[1]
  x2.span <- x2.range[2] - x2.range[1]
  y <- d$intensity
  y.range <- range(Re(y))

  # Normalizing data
  x1 <- (x1 - x1.range[1])/x1.span
  x2 <- (x2 - x2.range[1])/x2.span
  y <- cmplx2(rr = y$rr, ri = y$ri, ir = y$ir, ii = y$ii)
  y <- y/y.range[2]

  # Scaling and unpacking all parameters
  peaks <- peaks(object)
  n.peaks <- nrow(peaks)
  data.columns <- c('position', 'width', 'height', 'fraction.gauss')

  columns <- c("mixture", "species", "resonance")
  descriptors <- as.matrix(peaks[, which(colnames(peaks) %in% columns)])
  descriptors <- apply(descriptors, 1, paste, collapse = "-")

  
  dimensions <- peaks$dimension 

  peaks <- list(par = peaks[, data.columns], 
                lb = bounds(object)$lower$peaks[, data.columns], 
                ub = bounds(object)$upper$peaks[, data.columns])

  for (name in names(peaks)) {
    logic <- dimensions == "direct"

    peaks[[name]]$position <- 
      ifelse(logic, (peaks[[name]]$position - x1.range[1])/x1.span,
                    (peaks[[name]]$position - x2.range[1])/x2.span)

    peaks[[name]]$width <- 
      ifelse(logic, peaks[[name]]$width/direct.sf/x1.span,
                    peaks[[name]]$width/indirect.sf/x2.span)

    peaks[[name]]$height <- peaks[[name]]$height/sqrt(y.range[2])

    peaks[[name]] <- as.vector(t(as.matrix(peaks[[name]])))
  }

  #---------------------------------------
  # Scaling and unpacking baseline/phase terms
  
  # Scaling knots in line with x
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
  if ( n.phase == 3 ) {
    phase[1] <- phase[1] + phase[2]*x1.range[1] + phase[3]*x1.range[1]
    phase[2] <- (phase[2] - phase[1])*x1.span
    phase[3] <- (phase[3] - phase[1])*x2.span
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
      baseline[[name]]$rr/y.range[2], 
      baseline[[name]]$ri/y.range[2], 
      baseline[[name]]$ir/y.range[2], 
      baseline[[name]]$ii/y.range[2], 
      phase[[name]]
    )
  }

  #---------------------------------------
  # Generating constraint lists

  constraints <- parse_constraints(object, x1.span, x2.span)

  eq.constraints <- constraints[[1]]
  ineq.constraints <- constraints[[2]]

  #---------------------------------------
  # Performing the fit

  # Flattening constraints
  eq.constraints <- unlist(lapply(eq.constraints, function (x) c(x, NaN)))
  eq.constraints <- eq.constraints[-length(eq.constraints)]

  ineq.constraints <- unlist(lapply(ineq.constraints, function (x) c(x, NaN)))
  ineq.constraints <- ineq.constraints[-length(ineq.constraints)]

  # Encoding resonance/dimension information
  i.res <- as.integer(factor(descriptors))
  i.res <- rep(i.res - 1, each = 4)

  i.dim <- as.integer(factor(dimensions, levels = c('direct', 'indirect')))
  i.dim <- rep(i.dim - 1, each = 4)

  print(par)

  start.time <- proc.time()
  out <- .Call("fit_2d_wrapper", 
    x_direct = as.double(x1), 
    x_indirect = as.double(x2), 
    y = as.double(c(y$rr, y$ri, y$ir, y$ii)), 
    resonances = as.integer(i.res),
    dimensions = as.integer(i.dim),
    knots = as.double(knots),
    p = as.double(par$par), 
    lb = as.double(par$lb), 
    ub = as.double(par$ub), 
    n = as.integer(length(x1)),
    nl = as.integer(n.peaks*4),
    nb = as.integer(n.baseline),
    np = as.integer(n.phase),
    nk = as.integer(length(knots)),
    eq = as.double(eq.constraints),
    iq = as.double(ineq.constraints),
    neq = as.integer(length(eq.constraints)),
    niq = as.integer(length(ineq.constraints))
  )
  object@time <- as.numeric(proc.time() - start.time)[3]

  par$par <- out
  print(out)

  #---------------------------------------
  # Unpacking and rescaling parameters
  
  # Starting with peaks
  peaks <- peaks(object)
  new.peaks <- matrix(par$par[1:(n.peaks*4)], ncol = 4, byrow = TRUE)
  peaks[, data.columns] <- new.peaks

  peaks$position <- ifelse(
    dimensions == "direct",
    peaks$position*x1.span + x1.range[1],
    peaks$position*x2.span + x2.range[1])

  peaks$width <- ifelse(
    dimensions == "direct",
    peaks$width*direct.sf*x1.span,
    peaks$width*indirect.sf*x2.span)

  peaks$height <- peaks$height*sqrt(y.range[2])

  peaks(object) <- peaks

  # Then baseline
  index.rr <- (n.peaks*4 + 1):(n.peaks*4 + n.baseline)
  index.ri <- index.rr + n.baseline
  index.ir <- index.ri + n.baseline
  index.ii <- index.ir + n.baseline

  if ( n.baseline > 0 ) {
    object@baseline <- cmplx2(rr = par$par[index.rr]*y.range[2],
                              ri = par$par[index.ri]*y.range[2],
                              ir = par$par[index.ir]*y.range[2],
                              ii = par$par[index.ii]*y.range[2])
  }

  # And phase
  index <- (n.peaks*4 + n.baseline*4 + 1):(n.peaks*4 + n.baseline*4 + n.phase)
  if ( n.phase > 0 ) {
    phase <- par$par[index]

    # 1st order phase coefficients must be adapted back to the global scale
    if ( n.phase == 3 ) {
      phase.left <- phase[1]
      phase.right.1 <- phase[1] + phase[2]
      phase.right.2 <- phase[1] + phase[3]
      phase[2] <- (phase.right.1 - phase.left)/x1.span
      phase[3] <- (phase.right.2 - phase.left)/x2.span
      phase[1] <- phase.left - phase[2]*x1.range[1] - phase[3]*x2.range[1]
    }

    object@phase <- phase
  }

  object

})




#==============================================================================>
#  Display function
#==============================================================================>



#------------------------------------------------------------------------------
#' Display NMRFit2D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRFit2D", 
  function(object) {

    # Generating compiled data frames
    peaks <- peaks(object)
    baseline <- baseline(object)
    knots <- knots(object)
    phase <- phase(object)
    lower <- lower_bounds(object)
    upper <- upper_bounds(object)
    couplings <- couplings(object)

    cat('An object of NMRFit2D class\n\n')

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
setMethod("bounds", "NMRFit2D", 
  function(object) {
    # Extracting peak bounds from species using the NMRMixture2D signature
    bounds <- selectMethod("bounds", signature="NMRMixture2D")(object)

    # Baseline and phase
    lower.baseline <- object@lower.bounds$baseline
    if ( length(lower.baseline) == 0 ) {
      lower.baseline <- cmplx2(rr = -Inf, ri = -Inf, ir = -Inf, ii = -Inf)
    } 
    upper.baseline <- object@upper.bounds$baseline
    if ( length(upper.baseline) == 0 ) {
      upper.baseline <- cmplx2(rr = Inf, ri = Inf, ir = Inf, ii = Inf)
    }

    lower.phase <- object@lower.bounds$phase
    if ( length(lower.phase) == 0 ) {
      lower.phase <- -Inf
    }
    upper.phase <- object@upper.bounds$phase
    if ( length(upper.phase) == 0 ) {
      upper.phase <- Inf
    }

    # Outputting
    list(lower = list(peaks = bounds$lower, baseline = lower.baseline, 
                      phase = lower.phase), 
         upper = list(peaks = bounds$upper, baseline = upper.baseline,
                      phase = upper.phase))
  })



#------------------------------------------------------------------------------
# Baseline


#' @rdname baseline
#' @export
setMethod("baseline", "NMRFit2D", 
  function(object) object@baseline
)


#------------------------------------------------------------------------------
# Knots



#' @rdname knots
#' @export
setMethod("knots", "NMRFit2D", 
  function(object) object@knots
)


#------------------------------------------------------------------------------
# Phase

#' @rdname phase
#' @export
setMethod("phase", "NMRFit2D", 
  function(object) object@phase
)
