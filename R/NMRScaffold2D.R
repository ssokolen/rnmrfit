# Definition of a super-class for 2D resonance data.



#==============================================================================>
#  NMRScaffold2D -- super-class for Resonance, Species, Mixture and Fit
#==============================================================================>



#------------------------------------------------------------------------------
#' Super-class for all 2D peak descriptions.
#' 
#' This class is not meant to be used directly. Instead, it provides a common
#' framework for methods around visualization and inspection of NMRScaffold2D,
#' NMRSpecies2D, NMRMixture2D, and NMRFit2D objects. Essentially, all of the 2D
#' methods merely wrap around their 1D counterparts, meaning that the same
#' basic approach can be used for all of them.
#'
#' @slot direct An NMRScaffold1D object.
#' @slot indirect An NMRScaffold1D object.
#' 
#' @name NMRScaffold2D-class
#' @export
NMRScaffold2D <- setClass("NMRScaffold2D",
  contains = "VIRTUAL",
)



#==============================================================================>
# Access to direct and indirect components
#==============================================================================>



#------------------------------------------------------------------------------
# Direct

#---------------------------------------
#' Get 1D projection of direct dimension 
#' 
#' Generic convenience method to extract all direct dimension components.
#' 
#' @param object An NMRScaffold2D object.
#' 
#' @name direct
setGeneric("direct", 
  function(object) standardGeneric("direct")
)

#---------------------------------------
#' Set direct dimension resonance object
#' 
#' Generic convenience method to set the direct dimension resonance component
#' of an NMRScaffold2D object.
#' 
#' @param object An NMRScaffold2D object.
#' @param value A valid NMRScaffold1D object.
#' 
#' @name direct-set
#' @export
setGeneric("direct<-", 
  function(object, value) standardGeneric("direct<-")
)



#------------------------------------------------------------------------------
# Indirect

#---------------------------------------
#' Get 1D projection of indirect dimension 
#' 
#' Generic convenience method to extract all indirect dimension components.
#' 
#' @param object An NMRScaffold2D object.
#' 
#' @name indirect 
setGeneric("indirect", 
  function(object) standardGeneric("indirect")
)

#---------------------------------------
#' Set indirect dimension resonance object
#' 
#' Generic convenience method to set the indirect dimension resonance component
#' of an NMRScaffold2D object.
#' 
#' @param object An NMRScaffold2D object.
#' @param value A valid NMRScaffold1D object.
#' 
#' @name indirect-set
#' @export
setGeneric("indirect<-", 
  function(object, value) standardGeneric("indirect<-")
)



#==============================================================================>
#  Display function
#==============================================================================>



#------------------------------------------------------------------------------
#' Display NMRScaffold2D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRScaffold2D", 
  function(object) {

    direct <- direct(object)
    indirect <- indirect(object)
    id <- ifelse( 'id' %in% slotNames(object), sprintf('(%s)', id(object)), '')

    cat('\n\n')
    msg <- sprintf('\nAn object of %s class (%s)\n\n', class(object), id)
    cat(strrep('=', nchar(msg) - 3))
    cat(msg)

    cat('------------------------\n')
    cat('In the direct dimension:\n')
    cat('------------------------\n\n')
    show(direct)

    cat('\n\n')
    cat('--------------------------\n')
    cat('In the indirect dimension:\n')
    cat('--------------------------\n\n')
    show(indirect)

})
