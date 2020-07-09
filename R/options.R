# Definition of global options

#---------------------------------------
# The general options list

opts <- list(

  "sf" = list(
    .value = NULL,
    .length = 1,
    .class = "numeric"
  ),
  
  "baseline" = list(
    .value = list(order = 3, n.knots = 0),
    .length = 2,
    .class = "list",
    .validate = function (x) { 
      all(names(x) %in% c('degree', 'n.knots')) &&
      all(lapply(x, class) == 'numeric')
    }
  ),
  
  "phase" = list(
    .value = list(order = 0),
    .length = 1,
    .class = "list",
    .validate = function (x) { 
      ( names(x) == 'order' ) &&
      all(lapply(x, class) == 'numeric')
    } 
  ),

  "fit" = list(
    .value = list(
      opts = list(),
      init = function (object, ...) {
        object <- initialize_heights(object, nmrdata = object@nmrdata)
        object <- set_conservative_bounds(object, nmrdata = object@nmrdata)
        object
      }
    ),
    .length = 2,
    .class = "list",
    .validate = function (x) {
      all(names(x) %in% c('opts', 'init')) &&
      ( class(x[['init']]) == 'function' )
    }
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


