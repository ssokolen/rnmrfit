#------------------------------------------------------------------------------
check_general_bounds <- function(object, par, tolerance) {

  testthat_tolerance(tolerance)
  n <- nrow(peaks(object))

  for ( p in names(par) ) {
    if ( p %in% c("position", "width", "height", "fraction.gauss")) {

      expect_equal(lower_bounds(object)[[p]], rep(par[[p]][[1]], n))
      expect_equal(upper_bounds(object)[[p]], rep(par[[p]][[2]], n))
    }
  }
}

#------------------------------------------------------------------------------
check_offset_bounds <- function(object, par, tolerance) {

  testthat_tolerance(tolerance)
  n <- nrow(peaks(object))

  for ( p in names(par) ) {
    if ( p %in% c("position", "width", "height", "fraction.gauss")) {

      lower <- par[[p]][[1]]
      upper <- par[[p]][[2]]

      if ( par$relative ) {
        lower <- lower*peaks(object)[[p]]
        upper <- upper*peaks(object)[[p]]
      } 

      expect_equal(lower_bounds(object)[[p]], peaks(object)[[p]] + lower)
      expect_equal(upper_bounds(object)[[p]], peaks(object)[[p]] + upper)
    }
  }
}
