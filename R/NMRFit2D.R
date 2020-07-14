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

    peaks[[name]]$height <- peaks[[name]]$height/y.range[2]

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
  i.res <- as.integer(factor(descriptors)) - 1
  i.res <- rep(i.res, each = 4)

  i.dim <- as.integer(factor(dimensions, 
                             levels = c('direct', 'indirect'))) - 1
  i.dim <- rep(i.dim, each = 4)

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

  peaks$height <- peaks$height*y.range[2]

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
  index <- (n.peaks*4 + n.baseline*2 + 1):(n.peaks*4 + n.baseline*2 + n.phase)
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



#========================================================================>
#  Lineshape and area calculations
#========================================================================>



# TODO -- still from 1D

#' @rdname f_baseline
#' @export
setMethod("f_baseline", "NMRFit2D",
  function(object, components = 'r/i') {

    # Defining which components to return
    return.r <- grepl('r', tolower(components))
    return.i <- grepl('i', tolower(components))

    err <- '"components" must have at least one of either "r" or "i"'
    if ( return.r && return.i ) f_out <- function(y) {y}
    else if ( return.r ) f_out <- function(y) {Re(y)}
    else if ( return.i ) f_out <- function(y) {Im(y)}
    else stop(err)

    #knots <- object@knots
    #baseline <- object@baseline
    #order <- length(baseline) - length(knots) - 1

    # If there are no baseline parameters, return a dummy functions
    if ( TRUE ) { #length(baseline) == 0 ) {
      function(x1, x2) {
        zeros <- rep(0, length(x1))
        f_out(cmplx2(rr = zeros, ri = zeros, ir = zeros, ii = zeros))
      }
    }
    # Otherwise, generating actual baseline function
    else {
      function(x) {
        basis <- splines::bs(x, degree = order, knots = knots, intercept = TRUE)
        y <- cmplx1(r = as.vector(basis %*% matrix(Re(baseline), ncol = 1)),
                    i = as.vector(basis %*% matrix(Im(baseline), ncol = 1)))
        f_out(y)
      }
    }
  })



#==============================================================================>
# Plotting  
#==============================================================================>



#------------------------------------------------------------------------------
#' Plot NMRFit2D object
#' 
#' Generates an interactive plot object using the plotly package.
#' 
#' Convenience function that generates a graphical representation of the fit.
#' The original data is plotted as a black line, the fit is plotted in red, the
#' baseline is plotted in blue, the residual in red. The fit can be plotted as
#' a composite of all the peaks, or individually.
#' 
#' @param x An NMRFit2D object.
#' @param components One of 'rr', 'ii', 'ir', or 'ri' to specify various real
#'                   and imaginary components. 3D plot subplots are not
#'                   currently supported.
#' @param sum.level One of either 'all', 'species', 'resonance', to specify
#'                  whether all peaks should be summed together.
#' @param sum.baseline TRUE to add the baseline to each fit.
#' @param apply.phase TRUE to apply the calculated phase to the data.
#' 
#' @return A ggplot2 plot.
#' 
#' @export
plot.NMRFit2D <- function(x, components = 'rr', apply.phase = TRUE,  
                          sum.level = 'species', sum.baseline = TRUE) { 

  # Unlike 1D, 2D plots are limited to a single component
  components <- tolower(components)
  valid.components <- c('rr', 'ri', 'ir', 'ii')
  err <- '"components" must be one of "rr", "ri", "ir", or "ii".'
  if (! components %in% valid.components ) stop(err)

  # Selecting plot label
  if ( grepl('rr', components) ) {
    label <- "Real"
  } else if ( grepl('ri', components) )  {
    label <- "Real/Imaginary"
  } else if ( grepl('ir', components) ) {
    label <- "Imaginary/Real"
  } else {
    label <- "Imaginary"
  }

  #---------------------------------------
  # Calculating all required values

  nmrdata <- x@nmrdata 

  if ( apply.phase ) {
    nmrdata <- apply_phase(nmrdata, x@phase, degrees = FALSE)
  }

  # The original data
  d <- nmrdata@processed %>%
    arrange(direct.shift)
  direct.shift <- d$direct.shift 
  indirect.shift <- d$indirect.shift 

  # TODO
  y.data <- d$intensity$rr

  # The overall fit
  sf <- c(get_parameter(x@nmrdata, 'sfo1', 'acqus', dimension = "direct"),
          get_parameter(x@nmrdata, 'sfo1', 'acqus', dimension = "indirect"))

  f <- f_lineshape(x, sf = sf, sum.peaks = TRUE, components = components)
  y.fit <- f(direct.shift, indirect.shift)

  # The baseline
  f <- f_baseline(x, components = components)
  y.baseline <- f(direct.shift, indirect.shift)

  # The residual
  y.residual <- y.data - y.fit - y.baseline

  # All individual fits
  y.fit.all <- values(x, direct.shift, indirect.shift, sf = sf,
                      sum.peaks = FALSE, sum.baseline = FALSE,
                      components = components)

  # Generating grouped fits based on sum.level. The output is a list of
  # of data.frames with names that will be plotted one at a time

  # If everything is to be summed, generate frame from overall fit data
  if ( sum.level == 'all' ) {
    d <- data.frame(
      direct.shift = direct.shift, 
      indirect.shift = indirect.shift, 
      intensity = y.fit
    )
    frames <- list('Fit' = d)
  }
  else {
    err <- '"sum.level" must be one of "all", "species", or "resonance"'
    if ( sum.level == 'species' ) columns <- 'species'
    else if ( sum.level == 'resonance' ) columns <- c('species', 'resonance')
    else stop(err)
  
    # Tacking on chemical shift as a grouping column
    all.columns <- c(columns, 'direct.shift', 'indirect.shift')

    d <- y.fit.all %>%
      group_by_at(all.columns) %>%
      summarize(intensity = sum(intensity)) %>%
      ungroup()

    d$id <- apply(d[, columns], 1, paste, collapse = '-')

    d <- select(d, id, direct.shift, intensity)

    frames <- by(d, d$id, identity)
  }


  #---------------------------------------
  # Defining basic plot functions

  # Setting legend options
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

  # This function initializes the overall plot object by drawing a single
  # line with the colour and name of choice
  f_init <- function(z, color, name) {

    x <- direct.shift
    y <- indirect.shift
    groups <- unique(y)
    n <- length(groups)

    # Drawing separate lines for each indirect dimension
    index <- y == groups[1]
    p <- plot_ly(x = x[index], y = y[index], z = z[index],
                 name = I(name), color = I(color),
                 type = 'scatter3d', mode = 'lines', legendgroup = 1)

    # Looping over the rest
    for ( value in groups[-1] ) {
      index <- y == value
      p <- p %>%
        add_trace(x = x[index], y = y[index], z = z[index],
                  name = I(name), color = I(color),
                  type = 'scatter3d', mode = 'lines', legendgroup = 1,
                  showlegend = FALSE)
    }

    p %>%
      layout(scene = scene)
  }

  # This functions adds a new line to an existing plot object
  f_add <- function(p, z, color, name, group, showlegend = TRUE) {

    x <- direct.shift
    y <- indirect.shift
    groups <- unique(y)
    n <- length(groups)

    index <- y == groups[1]
    p <- p %>%
      add_trace(x = x[index], y = y[index], z = z[index],
                name = I(name), color = I(color),
                type = 'scatter3d', mode = 'lines', legendgroup = group)

    for ( value in groups[-1] ) {
      index <- y == value
      p <- p %>%
        add_trace(x = x[index], y = y[index], z = z[index],
                  name = I(name), color = I(color),
                  type = 'scatter3d', mode = 'lines', legendgroup = group,
                  showlegend = FALSE)
    }

    p
  }

  #---------------------------------------
  # Building up the plot elements

  # Initializing plots
  p <- f_init(y.data, 'black', label)

  # Adding baseline
  p <- f_add(p, y.baseline, 'blue', 'Baseline', 2) 

  # Adding residual
  p <- f_add(p, y.residual, 'green', 'Residual', 3) 

  # Looping through each previously defined peak grouping
  for ( i in 1:length(frames) ) {

    y <- frames[[i]]$intensity
    if ( sum.baseline ) y <- y + y.baseline

    id <- names(frames)[i] 

    p <- f_add(p, y, 'red', id, id) 
  }

  p
}

setMethod("plot", "NMRFit2D", plot.NMRFit2D)
