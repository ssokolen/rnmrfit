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
#' @name NMRFit1D-class
#' @export
NMRFit1D <- setClass("NMRFit1D",
  contains = 'NMRMixture1D',
  slots = c(
    nmrdata = 'NMRData1D',
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
#' NMRFit1D validity test
#'
validNMRFit1D <- function(object) {

  nmrdata <- object@nmrdata
  knots <- object@knots
  baseline <- object@baseline
  phase <- object@phase
  bounds <- object@bounds

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Checking nmrdata
  if ( (class(nmrdata) != 'NMRData1D') || (! validObject(nmrdata))  ) {

      valid <- FALSE
      new.err <- '"nmrdata" must be a valid NMRData1D object.'
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
setValidity("NMRFit1D", validNMRFit1D)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRFit1D object
#' 
#' Generates an NMRFit1D object based on a list of NMRSpecies1D objects or
#' other objects that can be converted to NMRSpecies1D objects. See
#' ?nmrresonance_1d and ?nmrspecies_1d for more details about this conversion.
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
#' @param species A list of NMRSpecies1D objects or other objects that can be
#'                converted to NMRSpecies1D objects. See ?nmrresonance_1d and
#'                ?nmrspecies_1d for more details about this conversion. If list
#'                elements are named, these names will be use to replace
#'                species ids.
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
#'                    polynomial. Only 0 and 1st order terms are typically used
#'                    in practice, but a higher order correction is possible.
#' @param delay.fit FALSE to immediately run least squares optimization after
#'                  the NMRFit1D object is initialized, TRUE to skip the
#'                  optimization, enabling more customization. The fit can be
#'                  run manually using fit().
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmroptions$direct$sf, but an override can be
#'           provided here.
#' @param init An initialization, function that takes an NMRFit1D object and
#'             returns a modified NMRFit1D object. Use the "identity" function
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
#' @param ... Options passed to nmrspecies_1d if conversion has to be performed.
#'            See ?nmrspecies_1d for more details.
#' 
#' @return An NMRFit1D object.
#' 
#' @export
nmrfit_1d <- function(
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
  mixture <- nmrmixture_1d(species, ...)

  #---------------------------------------
  # Checking nmrdata

  if ( (class(nmrdata) != 'NMRData1D') || (! validObject(nmrdata)) ) {
    err <- '"nmrdata" must be a valid NMRData1D object.'
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
  out <- new('NMRFit1D', mixture, nmrdata = nmrdata,
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



#------------------------------------------------------------------------------
#' Fit NMR data to peaks.
#' 
#' Given an NMRFit1D object, this function performs a least squares fit on the
#' data and updates the peak parameters. One of the most important things to
#' consider is that any peak found outside the data range of the NMRData1D
#' objects is automatically dropped from the resulting object. Exactly how this
#' happens can be modified using the exclusion.level and exclusion.notification
#' arguments. The exclusion.level parameter determines which part of the
#' overall species to exclude if any of its peaks fall outside the data range:
#' either 'species' for whole species, 'resonance' for just a subset of the
#' species and 'peak' to ignore resonance/species blocks and exclude by
#' specific peak alone. The exclusion.notification parameter determines how to
#' report when peaks are found to be outside the data range: either 'none' to
#' give no notice, 'message' to issue a message, 'warning' to issue a warning,
#' or 'stop' to issue an error.
#' 
#' @param object An NMRFit1D object.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmroptions$direct$sf = ..., but an override can be
#'           provided here.
#' @param init An initialization function that takes NMRFit1D object and returns
#'             a modified NMRFit1D object. Use the "identity" function to
#'             override the default initialization in the
#'             nmroptions$fit$init. Note that all arguments provided to the
#'             fit() function are also passed on the init() function.
#' @param opts A list of NLOPT fit options to override the default options in
#'             the nmroptionsd$fit$opts.
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
#' 
#' @name fit
#' @export
setGeneric("fit", 
  function(object, ...) {
    standardGeneric("fit")
})

#' @rdname fit
#' @export
setMethod("fit", "NMRFit1D",
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
      phase[2] <- phase[2]/x.span
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
              new.constraint <- c(3, -leeway, 
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
    fit_lineshape_1d(x, y, par$par, par$lb, par$ub, basis, 
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
      
      # 1st order phase coefficients must be adapted back to the global scale
      if ( n.phase == 2 ) {
        phase[2] <- phase[2]*x.span
        phase[1] <- phase[1] - phase[2]*x.range[1]
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
    bounds <- bounds(object)
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
setMethod("bounds", "NMRFit1D", 
  function(object) {
    # Extracting peak bounds from species using the NMRMixture1D signature
    bounds <- selectMethod("bounds", signature="NMRMixture1D")(object)

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

#---------------------------------------
#' Set object baseline
#' 
#' Generic convenience method to set the baseline definition of an NMRFit1D
#' object.
#' 
#' @param object An NMRFit1D object.
#' @param value A vector of spectra intensity values used to construct the
#'              baseline b-spline. The order of the baseline spline function is
#'              calculated automatically based on the length of the "knots" slot
#'              and the length of the baseline vector, where order =
#'              legnth(baseline) - length(knots).
#' 
#' @name baseline-set
#' @export
setGeneric("baseline<-", 
  function(object, value) standardGeneric("baseline<-")
)

#' @rdname baseline-set
#' @export
setReplaceMethod("baseline", "NMRFit1D",
  function(object, value) {

    if ( class(value) == 'numeric' ) {
      wrn <- paste('Applying the same baseline parameters to both real and',
                   'imaginary components.')
      warning(wrn)
      value <- complex(re = value, im = value)
    }

    object@baseline <- value
    validObject(object)
    object 
  })



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

#---------------------------------------
#' Set object baseline knots
#' 
#' Generic convenience method to set the baseline knots definition of an
#' NMRFit1D object.
#' 
#' @param object An NMRFit1D object.
#' @param value A vector of chemical shifts used to designate the internal
#'              baseline b-spline knots. Note that the number of internal knots
#'              impacts the length of the baseline and if the length of knots is
#'              changed, current baseline values may not feasible with the new
#'              knots.
#' 
#' @name knots-set
#' @export
setGeneric("knots<-", 
  function(object, value) standardGeneric("knots<-"))

#' @rdname knots-set
#' @export
setReplaceMethod("knots", "NMRFit1D",
  function(object, value) {

    # If the knot length changes, baseline parameters have to change
    if ( length(object@knots) != length(value) ) {
      
      # Generating y values from current baseline parameters
      x <- object@nmrdata@processed$direct.shift
      k1 <- object@knots
      b1 <- object@baseline
      n1 <- length(b1) - length(k1)
      X1 <- bs(x, degree = n1, knots = k1)
      y1 <- X1 %*% b1

      # Generating new baseline values from new basis
      k2 <- value
      X2 <- bs(x, degree = n1, knots = k2)
      b2 <- solve(t(X2) %*% X2) %*% t(X2) %*% y1
      y2 <- X2 %*% b2

      # If the resulting change in baseline parameters resulted in a change
      # to baseline values, issue a warning to that effect
      wrn <- paste('New knot values can not be used to represent current',
                   'baseline. New baseline parameters will be generated using',
                   'a least-squares fit.')
      if ( any( ( Re(y2-y1) > 1e-6 ) | ( Im(y2-y1) > 1e-6) ) ) warning(wrn)

      object@baseline <- as.vector(b2)
    }

    object@knots <- value
    validObject(object)
    object 
  })



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

#---------------------------------------
#' Set object phase
#' 
#' Generic convenience method to set the phase definition of an NMRFit1D
#' object.
#' 
#' @param object An NMRFit1D object.
#' @param value A vector of phase polynomial values in increasing order (in
#'              radians). So c(pi/4, pi/8) represents a 0 order phase correction
#'              of 45 degrees and a 1st order correction of 22.5 degrees
#'              (referenced to 0 ppm).
#' 
#' @name phase-set
#' @export
setGeneric("phase<-", 
  function(object, value) standardGeneric("phase<-"))

#' @rdname phase-set
#' @export
setReplaceMethod("phase", "NMRResonance1D",
  function(object, value) {
    object@phase <- value
    validObject(object)
    object 
  })



#==============================================================================>
#  Bounds
#==============================================================================>


#------------------------------------------------------------------------------
#' @rdname set_general_bounds
#' @export
setMethod("set_general_bounds", "NMRFit1D",
  function(object, ..., nmrdata = NULL, widen = FALSE,
           baseline = NULL, phase = NULL) {
  
  # First, propagating bounds to component species
  object@species <- lapply(object@species, set_general_bounds, ...,
                           nmrdata = nmrdata, widen = widen)

  #---------------------------------------
  # Then dealing with baseline and phase

  # Temporarily splitting real and imaginary baselines into two different values 
  if ( class(baseline) == 'numeric' ) {
    re.baseline <- baseline
    im.baseline <- NULL
  }
  else if ( class(baseline) == 'complex' ) {
    re.baseline <- Re(baseline)
    im.baseline <- Im(baseline)
  }
  else {
    re.baseline <- NULL
    im.baseline <- NULL
  }

  # Scaling baseline if nmrdata provided
  if (! is.null(nmrdata) ) {
    processed <- nmrdata@processed
    re.y.range <- max(Re(processed$intensity)) - min(Re(processed$intensity))
    im.y.range <- max(Im(processed$intensity)) - min(Im(processed$intensity))

    re.baseline <- re.baseline * re.y.range
    im.baseline <- im.baseline * im.y.range
  }

  # Defining a bound check function
  .check_bounds <- function(bounds) {

    if ( length(bounds) != 2 ) {
      err <- paste("All bounds must be vectors of two elements consisting",
                   "of a lower and upper bound.")
      stop(err)
    }

    if ( bounds[1] > bounds[2] ) {
      err <- paste("Lower bound must be smaller than upper bound.",
                   "Proceeding with current constraints will result in a",
                   "fit error.")
      warning(err)
    }

  }

  # Applying
  lower <- bounds(object)$lower
  upper <- bounds(object)$upper

  lower$re.baseline <- Re(lower$baseline)
  lower$im.baseline <- Im(lower$baseline)

  upper$re.baseline <- Re(upper$baseline)
  upper$im.baseline <- Im(upper$baseline)

  bounds = list(re.baseline = re.baseline, im.baseline = im.baseline,
                phase = phase)

  for ( parameter in names(bounds) ) {
    if ( length(bounds[[parameter]]) > 0 ) {
      .check_bounds(bounds[[parameter]])

      new <- bounds[[parameter]][1]
      if ( (new > lower[[parameter]]) || widen ) lower[[parameter]] <- new

      new <- bounds[[parameter]][2]
      if ( (new < upper[[parameter]]) || widen ) upper[[parameter]] <- new
    }
  }

  lower$baseline <- complex(re = lower$re.baseline, im = lower$im.baseline)
  upper$baseline <- complex(re = upper$re.baseline, im = upper$im.baseline)

  object@bounds$lower <- lower[c('baseline', 'phase')]
  object@bounds$upper <- upper[c('baseline', 'phase')]

  validObject(object)
  object
  })



#------------------------------------------------------------------------------
#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRFit1D",
  function(object, ...) {
    object@species <- lapply(object@species, set_offset_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRFit1D",
  function(object,  ..., nmrdata = NULL, widen = FALSE,
           baseline = TRUE, phase = TRUE) {
  
  # First, propagating bounds to component species
  object@species <- lapply(object@species, set_conservative_bounds, ...,
                           nmrdata = nmrdata, widen = widen)

  #---------------------------------------
  # Then dealing with baseline and phase

  # Adding baseline constraint if nmrdata is provided
  if ( (! is.null(nmrdata)) && baseline ) {
    object <- set_general_bounds(object, baseline = c(-0.5, 0.5),
                                 nmrdata = nmrdata, widen = widen)
  }

  # Phase constrain is applied regardless of nmrdata
  if ( phase ) {
    object <- set_general_bounds(object, phase = c(-pi/2, pi/2), widen = widen)
  }

  object
})



#------------------------------------------------------------------------------
#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRFit1D",
  function(object, ...) {
    object@species <- lapply(object@species, set_peak_type, ...)
    object
  })



#========================================================================>
#  Lineshape and area calculations
#========================================================================>



#------------------------------------------------------------------------
#' Generate baseline function
#' 
#' This is primarily an internal method that outputs a function that outputs
#' spectral intensity data of the fit baseline given a vector input of chemical
#' shifts.
#' 
#' @param object An NMRFit1D object.
#' @param components 'r/i' to output both real and imaginary data, 'r' to output
#'                   only real and 'i' to output only imaginary.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @return A function that outputs spectral intensity data of the fit baseline
#'         given a vector input of chemical shifts.
#' 
#' @name f_baseline
#' @export
setGeneric("f_baseline", 
  function(object, components = 'r/i', ...) standardGeneric("f_baseline")
)

#' @rdname f_baseline
#' @export
setMethod("f_baseline", "NMRFit1D",
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
        f_out(complex(re = zeros, im = zeros))
      }
    }
    # Otherwise, generating actual baseline function
    else {
      function(x) {
        basis <- splines::bs(x, degree = order, knots = knots, intercept = TRUE)
        y <- as.vector(basis %*% matrix(baseline, ncol = 1))
        f_out(y)
      }
    }
  })



#==============================================================================>
# Plotting  
#==============================================================================>



#------------------------------------------------------------------------------
#' Plot NMRFit1D object
#' 
#' Generates an interactive plot object using the plotly package.
#' 
#' Convenience function that generates a graphical representation of the fit.
#' The original data is plotted as a black line, the fit is plotted in red, the
#' baseline is plotted in blue, the residual in red. The fit can be plotted as
#' a composite of all the peaks, or individually.
#' 
#' @param x An NMRFit1D object.
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
plot.NMRFit1D <- function(x, components = 'r', apply.phase = TRUE,  
                          sum.level = 'species', sum.baseline = TRUE) { 

  #---------------------------------------
  # Calculating all required values

  # The original data
  d <- x@nmrdata@processed
  direct.shift <- d$direct.shift
  y.data <- d$intensity

  if ( apply.phase ) {
    y.data <- phase_spectrum(y.data, x@phase, degrees = FALSE)
  }

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
