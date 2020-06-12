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
    baseline = 'complex',
    phase = 'numeric',
    bounds = 'list',
    time = 'numeric'
  ),
  prototype = prototype(
    knots = numeric(0), 
    baseline = complex(re = rep(0, 4), im = rep(0, 4)),
    phase = c(0),
    bounds = list(lower = NULL, upper = NULL),
    time = c(0)
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRFit2D validity test
#'
validNMRFit2D <- function(object) {

  nmrdata <- object@nmrdata
  knots <- object@knots
  baseline <- object@baseline
  phase <- object@phase
  bounds <- object@bounds

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Checking nmrdata
  if ( (class(nmrdata) != 'NMRData2D') || (! validObject(nmrdata))  ) {

      valid <- FALSE
      new.err <- '"nmrdata" must be a valid NMRData2D object.'
      err <- c(err, new.err)

  }

  #---------------------------------------
  # Checking baseline length 
  if ( (length(baseline) > 0) && (length(baseline) <= (length(knots)+1)) ) {

      valid <- FALSE
      new.err <- paste('"baseline" vector length must be greater than the',
                       '"knots" vector length by at least 2 elements.')
      err <- c(err, new.err)

  }

  #---------------------------------------
  # Checking phase length 
  if ( length(phase) > 2  ) {

      valid <- FALSE
      new.err <- 'The phase correction term must be of length 2 or smaller.'
      err <- c(err, new.err)

  }



  #---------------------------------------
  # Checking that lower bounds match slots 
  valid.bounds <- c('baseline', 'phase')
  if (! is.null(bounds$lower) ) {

    logic <- identical(names(bounds$lower), valid.bounds)
    if (! logic ) {
      valid <- FALSE
      new.err <- sprintf('"bounds$lower" must have the following elements: %s',
                         paste(valid.bounds, collapse = ', '))
      err <- c(err, new.err)
    }

    logic.1 <- length(bounds$lower$baseline) %in% c(0, 1)
    logic.2 <- length(bounds$lower$phase) %in% c(0, 1) 
    if (! (logic.1 && logic.2) ) {
      valid <- FALSE
      new.err <- paste('"bounds$lower$baseline" and "bounds$lower$phase must',
                       'have a length of either zero or one, representing an',
                       'overall bound on baseline or phase correction.')
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking that upper bounds match slots
  if (! is.null(bounds$upper) ) {

    logic <- identical(names(bounds$upper), valid.bounds)
    if (! logic ) {
      valid <- FALSE
      new.err <- sprintf('"bounds$upper" must have the following elements: %s',
                         paste(valid.bounds, collapse = ', '))
      err <- c(err, new.err)
    }

    logic.1 <- length(bounds$upper$baseline) %in% c(0, 1)
    logic.2 <- length(bounds$upper$phase) %in% c(0, 1) 
    if (! (logic.1 && logic.2) ) {
      valid <- FALSE
      new.err <- paste('"bounds$upper$baseline" and "bounds$upper$phase must',
                       'have a length of either zero or one, representing an',
                       'overall bound on baseline or phase correction.')
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking the knots are all inside the boundaries
  direct.shift <- range(nmrdata@processed$direct.shift)
  if ( any((knots < direct.shift[1]) | (knots > direct.shift[2])) ) {

      wrn <- paste('It is recommended to keep "knots" values inside the',
                   'chemical shift range of the data.')
      warning(wrn)

  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRFit2D", validNMRFit2D)



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
#'                    polynomial. Only 0 and 1st order terms are typically used
#'                    in practice, but a higher order correction is possible.
#' @param delay.fit FALSE to immediately run least squares optimization after
#'                  the NMRFit2D object is initialized, TRUE to skip the
#'                  optimization, enabling more customization. The fit can be
#'                  run manually using fit().
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmroptions$direct$sf, but an override can be
#'           provided here.
#' @param init An initialization, function that takes an NMRFit2D object and
#'             returns a modified NMRFit2D object. Use the "identity" function
#'             to override the default initialization in the
#'             nmroptions$fit$init. Note that all arguments provided to the
#'             fit() function are also passed on the init() function.
#' @param opts A list of NLOPT fit options to override the default options in
#'             the nmroptions$fit$opts.
#' @param exclusion.level A string specifying what to do when peaks are found to
#'                        fall outside of the data range: either 'species' to
#'                        exclude the whole species to which the offending peak
#'                        belongs, 'resonance' to exclude the resonance to which
#'                        the offending peak belongs, or 'peak' to exclude just
#'                        the peak itself.
#' @param exclusion.notification A function specifying how to report when peaks
#'                               are found to be outside the data range: 'none'
#'                               to ignore, 'message' to issue a message,
#'                               'warning' to issue a warning, or 'stop' to
#'                               issue an error.
#' @param ... Options passed to nmrspecies_2d if conversion has to be performed.
#'            See ?nmrspecies_2d for more details.
#' 
#' @return An NMRFit2D object.
#' 
#' @export
nmrfit_2d <- function(
  species, nmrdata, baseline.order = nmroptions$direct$baseline$order,
  n.knots = nmroptions$direct$baseline$n.knots, 
  phase.order = nmroptions$direct$phase$order, 
  delay.fit = FALSE, sf = nmroptions$direct$sf,
  init = nmroptions$fit$init, opts = nmroptions$fit$opts, 
  exclusion.level = nmroptions$exclusion$level, 
  exclusion.notification = nmroptions$exclusion$notification, ...) {

  #---------------------------------------
  # Generating list of species 

  # Generate an NMRMixture object that will be used as a template for the fit
  mixture <- nmrmixture_2d(species, ...)

  #---------------------------------------
  # Checking nmrdata

  if ( (class(nmrdata) != 'NMRData2D') || (! validObject(nmrdata)) ) {
    err <- '"nmrdata" must be a valid NMRData2D object.'
    stop(err)
  }

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
    direct.shift <- range(nmrdata@processed$direct.shift)
    knots <- seq(direct.shift[1], direct.shift[2], length.out = n.knots)
  }

  # The initial value for the phase is always 0
  if ( phase.order == -1 ) phase <- numeric(0)
  else phase <- rep(0, phase.order + 1)

  # Initializing bounds
  bounds <- list(lower = list(baseline = complex(0), phase = numeric(0)),
                 upper = list(baseline = complex(0), phase = numeric(0)))

  #---------------------------------------
  # Resulting fit object
  out <- new('NMRFit2D', mixture, nmrdata = nmrdata,
                         knots = knots, baseline = baseline, phase = phase,
                         bounds = bounds)

  # If the fit is delayed, then return current object, otherwise run fit first
  if ( delay.fit ) out
  else fit(out, sf = sf, init = init, opts = opts, 
           exclusion.level = exclusion.level, 
           exclusion.notification = exclusion.notification)

}



#==============================================================================>
#  The main fit function (wrapping Rcpp code)
#==============================================================================>



#' @rdname fit
#' @export
setMethod("fit", "NMRFit2D",
  function(object, sf = nmroptions$direct$sf, init = nmroptions$fit$init, 
           opts = nmroptions$fit$opts, 
           exclusion.level = nmroptions$exclusion$level,
           exclusion.notification = nmroptions$exclusion$notification) {

    # First, run the initialization
    object <- init(object, sf = sf, init = init, opts = opts, 
                   exclusion.level = exclusion.level, 
                   exclusion.notification = exclusion.notification)

    # Unpacking some of the NMR data
    d <- processed(object@nmrdata)
    x <- d$direct.shift
    x.range <- range(x)
    x.span <- x.range[2] - x.range[1]
    y <- d$intensity
    y.range <- range(Re(y))

    # Normalizing data
    x <- (x - x.range[1])/x.span
    y <- complex(re = Re(y), im = Im(y))
    y <- y/y.range[2]

    # Baseline knots have to be scaled along with 

    # Exclude peaks in advance by tying into the update_peaks functions
    # (to consider full resonance/fit exclusion)
    peaks <- peaks(object)
    logic <- (peaks$position > x.range[1]) & (peaks$position < x.range[2])
    peaks <- peaks[logic, ]

    object2 <- update_peaks(object, peaks, exclusion.level = exclusion.level,
                            exclusion.notification = "none")

    peaks <- peaks(object2)
    n.peaks <- nrow(peaks)

    # Scaling and unpacking all parameters
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
    
    # Scaling knots in line with x
    knots <- knots(object)
    knots <- (knots - x.range[1])/x.span
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
      par[[name]] <- c(peaks[[name]], Re(baseline[[name]])/y.range[2], 
                        Im(baseline[[name]])/y.range[2], phase[[name]])
    }

    #---------------------------------------
    # Generating constraint lists
    # Each element in the list is an encoded constraint in the form of a vector.
    # The elements are decoded as follows:
    # 1: Code of either 0, 1, 2, 3, 4. 0 means position, 1 width, 2 height, 
    #    3 fraction.gauss, 4 area
    # 2: The equality or inequality target value. 
    # 3+: Integers corresponding to peak indexes in the overall parameter 
    #     vector. Positive vs negative values are treated differently depending
    #     on code. For positions, parameters with negative indexes are  
    #     subtracted while for width and area, all parameters with negative
    #     indexes are added up before dividing the positive ones.
    eq.constraints <- list()
    ineq.constraints <- list()

    # Redefining peaks from modified object
    peaks <- peaks(object)

    # To ensure that individual leeway values are considered, looping through
    # each species/resonance one at a time
    for (specie in object@species) {

      # At the species level, constraints are based on overall area sums
      leeway <- abs(specie@connections.leeway)

      # First dealing with conenctions between resonances
      connections <- specie@connections
      if ( nrow(connections) > 0 ) {
        
        for ( i in 1:nrow(connections) ) {
          resonance.1 <- connections$resonance.1[i]
          resonance.2 <- connections$resonance.2[i]
          ratio <- connections$area.ratio[i]

          logic <-(specie@id == peaks$species)
          logic.1 <- logic & (resonance.1 == peaks$resonance)
          logic.2 <- logic & (resonance.2 == peaks$resonance)

          if ( leeway == 0 ) {
            new.constraint <- c(4, ratio, which(logic.2), -which(logic.1))
            eq.constraints <- c(eq.constraints, list(new.constraint))
          }
          else {
            # The upper bounds
            new.constraint <- c(4, ratio*(1+leeway), 
                                which(logic.2), -which(logic.1))
            ineq.constraints <- c(ineq.constraints, list(new.constraint))

            # The lower bound
            new.constraint <- c(4, 1/(ratio*(1-leeway)), 
                                -which(logic.2), which(logic.1))
            ineq.constraints <- c(ineq.constraints, list(new.constraint))
          }
        }
      }

      # Then looping through each specific resonance within each species
      for (resonance in specie@resonances) {

        # At the species level, constraints are based on overall area sums
        id <- resonance@id
        position.leeway <- abs(resonance@couplings.leeway$position)
        width.leeway <- abs(resonance@couplings.leeway$width)
        fraction.leeway <- abs(resonance@couplings.leeway$fraction.gauss)
        area.leeway <- abs(resonance@couplings.leeway$area)

        # Only continue if there are couplings defined
        couplings <- resonance@couplings
        if ( nrow(couplings) > 0 ) {
          
          for ( i in 1:nrow(couplings) ) {
            peak.1 <- couplings$peak.1[i]
            peak.2 <- couplings$peak.2[i]
            diff <- couplings$position.difference[i]/x.span
            ratio <- couplings$area.ratio[i]

            logic <- (specie@id == peaks$species) & 
                       (resonance@id == peaks$resonance)
            logic.1 <- logic & (peak.1 == peaks$peak)
            logic.2 <- logic & (peak.2 == peaks$peak)

            # Each coupling constraint includes position, width, and area

            # First, the position
            leeway <- position.leeway
            if ( leeway == 0 ) {
              new.constraint <- c(0, diff, which(logic.2), -which(logic.1))
              eq.constraints <- c(eq.constraints, list(new.constraint))
            }
            else {
              # The upper bounds
              new.constraint <- c(0, diff*(1+leeway), 
                                  which(logic.2), -which(logic.1))
              ineq.constraints <- c(ineq.constraints, list(new.constraint))

              # The lower bound
              new.constraint <- c(0, -diff*(1-leeway), 
                                  -which(logic.2), which(logic.1))
              ineq.constraints <- c(ineq.constraints, list(new.constraint))
            }

            # Then, the width (same by default)
            leeway <- width.leeway
            if ( leeway == 0 ) {
              new.constraint <- c(1, 1, which(logic.2), -which(logic.1))
              eq.constraints <- c(eq.constraints, list(new.constraint))
            }
            else {
              # The upper bounds
              new.constraint <- c(1, (1+leeway), 
                                  which(logic.2), -which(logic.1))
              ineq.constraints <- c(ineq.constraints, list(new.constraint))

              # The lower bound
              new.constraint <- c(1, 1/(1-leeway), 
                                  -which(logic.2), which(logic.1))
              ineq.constraints <- c(ineq.constraints, list(new.constraint))
            }

            # Then, the fraction Gauss (same by default)
            leeway <- fraction.leeway
            if ( leeway == 0 ) {
              new.constraint <- c(3, 0, which(logic.2), -which(logic.1))
              eq.constraints <- c(eq.constraints, list(new.constraint))
            }
            else {
              # The upper bounds
              new.constraint <- c(3, leeway, 
                                  which(logic.2), -which(logic.1))
              ineq.constraints <- c(ineq.constraints, list(new.constraint))

              # The lower bound
              new.constraint <- c(3, leeway, 
                                  -which(logic.2), which(logic.1))
              ineq.constraints <- c(ineq.constraints, list(new.constraint))
            }

            # Finally, the area
            leeway <- area.leeway
            if ( leeway == 0 ) {
              new.constraint <- c(4, ratio, which(logic.2), -which(logic.1))
              eq.constraints <- c(eq.constraints, list(new.constraint))
            }
            else {
              # Inequality constraints around percentages can be problematic
              # since 10% of a 2x peak area difference doesn't directly map
              # to a 10% difference around the inverse 0.5x. As such, all
              # ratio are treated as being greater than 1, with indexes
              # flipped to accommodate.
              if ( ratio < 1 ) {
                ratio <- 1/ratio
                temp.logic <- logic.1
                logic.1 <- logic.2
                logic.2 <- temp.logic
              }

              # The upper bounds
              new.constraint <- c(4, ratio*(1+leeway), 
                                  which(logic.2), -which(logic.1))
              ineq.constraints <- c(ineq.constraints, list(new.constraint))

              # The lower bound
              new.constraint <- c(4, 1/(ratio*(1-leeway)), 
                                  -which(logic.2), which(logic.1))
              ineq.constraints <- c(ineq.constraints, list(new.constraint))
            }
          }
        }
      }
    }

    #---------------------------------------
    # Performing the fit

    # Generating baseline basis
    order <- n.baseline - length(knots) - 1
    if ( n.baseline > 0 ) {
      basis <- splines::bs(x, degree = order, knots = knots, intercept = TRUE)
    } else {
      # Creating a mock basis that will be ignored in the Rcpp code
      basis <- matrix(0, nrow = length(x), ncol = 1)
    }

    print(head(x))
    print(head(y))
    print(par)
    print(head(basis))
    print(eq.constraints)
    print(ineq.constraints)

    start.time <- proc.time()
    fit_lineshape_2d(x, y, par$par, par$lb, par$ub, basis, 
                     eq.constraints, ineq.constraints,
                     n.peaks, n.baseline, n.phase)
    object@time <- as.numeric(proc.time() - start.time)[3]

    #---------------------------------------
    # Unpacking and rescaling parameters
    
    # Starting with peaks
    peaks <- peaks(object)
    new.peaks <- matrix(par$par[1:(n.peaks*4)], ncol = 4, byrow = TRUE)
    peaks[, data.columns] <- new.peaks

    peaks$position <- peaks$position*x.span + x.range[1]
    peaks$width <- peaks$width*sf*x.span
    peaks$height <- peaks$height*y.range[2]

    object <- update_peaks(object, peaks, exclusion.level = exclusion.level,
                           exclusion.notification = exclusion.notification)

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

      print('Phase!')
      print(phase)
      print(x.range)
      
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
    bounds <- bounds(object)
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

    cat('Imaginary baseline: ')
    if ( length(baseline) > 0) {cat('\n'); print(round(Im(baseline), 4))}
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

    lower <- unlist(bounds$lower$peaks[ , columns])
    upper <- unlist(bounds$upper$peaks[ , columns])
    
    range <- paste('(', lower, ', ', upper, ')', sep = '')
    peaks[ , columns] <- range

    cat('Bounds (lower, upper):\n\n')
    print(peaks)
    cat('\n')   

    # Baseline bounds
    cat('Baseline and phase correction bounds:\n\n')

    cat('Real baseline: ')
    if ( length(bounds$lower$baseline) > 0) {
      lower <- round(Re(bounds$lower$baseline), 4)
      upper <- round(Re(bounds$upper$baseline), 4)
      
      range <- paste('(', lower, ', ', upper, ')\n', sep = '')
      cat(range)
    }
    else cat('None\n')

    cat('Imaginary baseline: ')
    if ( length(bounds$lower$baseline) > 0) {
      lower <- round(Im(bounds$lower$baseline), 4)
      upper <- round(Im(bounds$upper$baseline), 4)
      
      range <- paste('(', lower, ', ', upper, ')\n', sep = '')
      cat(range)
    }
    else cat('None\n')

    # Phase bounds
    cat('Phase (radians): ')

    if ( length(bounds$lower$phase) > 0)  {
      lower <- round(bounds$lower$phase, 4)
      upper <- round(bounds$upper$phase, 4)
      
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
    lower.baseline <- object@bounds$lower$baseline
    lower.baseline <- ifelse(length(lower.baseline) == 0, 
                             complex(re = -Inf, im = -Inf), lower.baseline)
    upper.baseline <- object@bounds$upper$baseline
    upper.baseline <- ifelse(length(upper.baseline) == 0, 
                             complex(re = Inf, im = Inf), upper.baseline)   

    lower.phase <- object@bounds$lower$phase
    lower.phase <- ifelse(length(lower.phase) == 0, -Inf, lower.phase)   
    upper.phase <- object@bounds$upper$phase
    upper.phase <- ifelse(length(upper.phase) == 0, Inf, upper.phase)   

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

    knots <- object@knots
    baseline <- object@baseline
    order <- length(baseline) - length(knots) - 1

    # If there are no baseline parameters, return a dummy functions
    if ( length(baseline) == 0 ) {
      function(x) {
        zeros <- rep(0, length(x))
        f_out(cmplx1(r = zeros, i = zeros))
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
#' @param components One of either 'r', 'i', or 'r/i' to include real, imaginary
#'                   or both components. If both components are selected, they
#'                   are displayed in separate subplots.
#' @param sum.level One of either 'all', 'species', 'resonance', 'peak' to
#'                  specify whether all peaks should be summed together the
#'                  peaks should be summed at a lower level.
#' @param sum.baseline TRUE to add the baseline to each fit.
#' @param apply.phase TRUE to apply the calculated phase to the data.
#' 
#' @return A ggplot2 plot.
#' 
#' @export
plot.NMRFit2D <- function(x, components = 'r', apply.phase = TRUE,  
                          sum.level = 'species', sum.baseline = TRUE) { 

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
  y.data <- d$intensity

  # The overall fit
  sf <- get_parameter(x@nmrdata, 'sfo1', 'acqus')
  if ( is.null(sf) ) sf <- nmroptions$direct$sf

  f <- f_lineshape(x, sf = sf, sum.peaks = TRUE)
  y.fit <- f(direct.shift)

  # The baseline
  f <- f_baseline(x)
  y.baseline <- f(direct.shift)

  # The residual
  y.residual <- y.data - y.fit - y.baseline

  # All individual fits
  y.fit.all <- values(x, direct.shift, sf = sf,
                      sum.peaks = FALSE, sum.baseline = FALSE)

  # Generating grouped fits based on sum.level. The output is a list of
  # of data.frames with names that will be plotted one at a time

  # If everything is to be summed, generate frame from overall fit data
  if ( sum.level == 'all' ) {
    d <- data.frame(direct.shift = direct.shift, intensity = y.fit)
    frames <- list('Fit' = d)
  }
  else {
    err <- '"sum.level" must be one of "all", "species", "resonance", or "peak"'
    if ( sum.level == 'species' ) columns <- 'species'
    else if ( sum.level == 'resonance' ) columns <- c('species', 'resonance')
    else if ( sum.level == 'peak' ) columns <- c('species', 'resonance', 'peak')
    else stop(err)
  
    # Tacking on direct.shift as a grouping column
    all.columns <- c(columns, 'direct.shift')

    # Converting cmplx1 to complex() for dplyr support
    y.fit.all$intensity <- vec_cast(y.fit.all$intensity, complex())
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

  # Note that the x values for each of the following functions is already
  # set as the direct shift of the data

  # This function initializes the overall plot object by drawing a single
  # line with the colour and name of choice
  f_init <- function(y, color, name) {
    p <- plot_ly(x = direct.shift, y = y, color = I(color), 
                 name = I(name), type = 'scatter', mode = 'lines',
                 legendgroup = 1) %>%
         layout(legend = legend.opts,
                xaxis = list(autorange = "reversed"))
  }

  # This functions adds a new line to an existing plot object
  f_add <- function(p, y, color, name, group, showlegend = TRUE) {
    p %>% 
      add_trace(x = direct.shift, y = y, color = I(color),
                name = I(name), type = 'scatter', mode = 'lines',
                legendgroup = group, showlegend = showlegend)
  }

  #---------------------------------------
  # Building up the plot elements

  # Initializing the plot list
  plots <- list()

  # Checking which components to plot
  re <- grepl('r', components)
  im <- grepl('i', components)

  # Initializing plots
  if ( re ) plots$r <- f_init(Re(y.data), 'black', 'Real')
  if ( im ) plots$i <- f_init(Im(y.data), 'grey', 'Imaginary')

  # Adding baseline
  if ( re ) plots$r <- f_add(plots$r, Re(y.baseline), 'blue', 'Baseline', 2) 
  if ( im ) plots$i <- f_add(plots$i, Im(y.baseline), 'blue', 'Baseline', 2)

  # Adding residual
  if ( re ) plots$r <- f_add(plots$r, Re(y.residual), 'green', 'Residual', 3) 
  if ( im ) plots$i <- f_add(plots$i, Re(y.residual), 'green', 'Residual', 3)

  # Looping through each previously defined peak grouping
  for ( i in 1:length(frames) ) {

    y <- frames[[i]]$intensity
    if ( sum.baseline ) y <- y + y.baseline

    id <- names(frames)[i] 

    if ( re ) plots$r <- f_add(plots$r, Re(y), 'red', id, id) 
    if ( im ) plots$i <- f_add(plots$i, Im(y), 'red', id, id)
  }

  if ( length(plots) == 0 ) NULL
  else subplot(plots, shareX = TRUE, shareY = TRUE, 
               nrows = min(length(plots), 2))
}
