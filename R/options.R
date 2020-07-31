# Definition of global options

#---------------------------------------
# The general options list

opts <- list(

  "sf" = list(
    .value = NULL,
    .length = 1,
    .class = "numeric"
  ),
  
  "baseline" = set_opt(
    order = list(
      .value = 3,
      .length = 1,
      .class = "numeric",
      .validate = function (x) { x %in% c(-1, 0:20) },
      .failed_msg = "Value must be -1 (disabled) or 0-20."
    ),
    n.knots = list(
      .value = 0,
      .length = 1,
      .class = "numeric",
      .validate = function (x) { x >= 0 },
      .failed_msg = "Value must be a positive integer."
    )
  ),
  
  "phase" = set_opt(
    order = list(
      .value = 0,
      .length = 1,
      .class = "numeric",
      .validate = function (x) { x %in% c(-1, 0, 1) },
      .failed_msg = "Value must be -1 (disabled), 0, or 1."
    )
  ),

  "fit" = set_opt(
    opts = set_opt(
      algorithm = list(
        .value = "SLSQP",
        .length = 1,
        .class = "character",
        .validate = function (x) { x %in% c("SLSQP", "ISRES") },
        .failed_msg = 'Only "SLSQP" and "ISRES" currently supported.'
      ),
      padding = list(
        .value = 0,
        .length = 1,
        .class = "numeric",
        .validate = function(x) {x >= 0},
        .failed_msg = "Value must be positive."
      ),
      xtol_rel = list(
        .value = 1e-4,
        .length = 1,
        .class = "numeric",
        .validate = function(x) {(x > 0) & (x < 1)},
        .failed_msg = "Value should be between 0 and 1."
      ),
      maxtime = list(
        .value = 0,
        .length = 1,
        .class = "numeric",
        .validate = function(x) {x >= 0},
        .failed_msg = "Value must be positive."
      )
    ),
    init = list(
      .value = function (object, ...) {
        object <- initialize_heights(object, nmrdata = object@nmrdata)
        object <- set_conservative_bounds(object, nmrdata = object@nmrdata)
        object
      },
      .class = "function"
    )
  )
)



#---------------------------------------
# Initializing global options

#' @export
nmroptions <- do.call(set_opt, opts)



#---------------------------------------
# Introducing direct and indirect sublists for sf
# (which may be expanded later as needed)

named.opts <- c('sf')
child_opts <- opts[named.opts]

for ( entry in named.opts ) {
  f_gen <- function (x) { 
    force(x)
    function () { nmroptions(x) } 
  }
  f <- f_gen(entry)
  child_opts[[entry]][['.value']] <- f
}

nmroptions_direct <- do.call(set_opt, child_opts)
nmroptions_indirect <- do.call(set_opt, child_opts)


nmroptions(direct = nmroptions_direct, ADD = TRUE)
nmroptions(indirect = nmroptions_indirect, ADD = TRUE)



#---------------------------------------
# Display sublist

nmroptions_display <- set_opt(
  digits = list(
    .value = 4,
    .length = 1,
    .class = "numeric"
  )
)

nmroptions(display = nmroptions_display, ADD = TRUE)


