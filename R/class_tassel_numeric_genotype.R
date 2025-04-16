## ----
#' @title
#' TasselNumericGenotype Class Definition
#' 
#' @description
#' Defines the `TasselNumericGenotype` class, which extends the 
#' `TasselGenotype` class. This class is used to represent numeric 
#' genotype data in the TASSEL 5 framework.
#' 
#' @name TasselNumericGenotype-class
#' @rdname TasselNumericGenotype-class
#' @exportClass TasselGenotype
setClass(
    Class = "TasselNumericGenotype",
    contains = "TasselGenotype"
)



# /// Methods (show) ////////////////////////////////////////////////

## ----
#' @title
#' Display Information for TasselNumericGenotype Object
#' 
#' @description
#' This method is used to display information about a 
#' \code{TasselNumericGenotype} object. It prints a summary of the 
#' genotype data, including the number of taxa, number of sites, and 
#' memory address of the Java object.
#' 
#' @details
#' The function \code{printNumGtDisp} is called internally to format 
#' and display the genotype data. The number of taxa and sites are 
#' retrieved from the Java reference object associated with the 
#' \code{TasselNumericGenotype} instance.
#'
#' @param object An object of class \code{TasselNumericGenotype}.
#'
#' @docType methods
#' @rdname TasselNumericGenotype-class
#' @aliases show,TasselNumericGenotype-method
setMethod("show", "TasselNumericGenotype", function(object) {
    printNumGtDisp(
        fgs    = object@dispData,
        nTaxa  = object@jRefObj$numberOfTaxa(),
        nSites = object@jRefObj$numberOfSites(),
        jMem   = object@jMemAddress
    )
})



# /// Methods (general) /////////////////////////////////////////////

## ----
#' @rdname javaRefObj
#' @export
setMethod(
    f = "javaRefObj",
    signature = signature(object = "TasselNumericGenotype"),
    definition = function(object) {
        return(object@jRefObj)
    }
)


