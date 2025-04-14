## ----
#' @name TasselNumericGenotype-class
#' @rdname TasselNumericGenotype-class
#' @exportClass TasselGenotype
setClass(
    Class = "TasselNumericGenotype",
    contains = "TasselGenotype"
)



# /// Methods (show) ////////////////////////////////////////////////

## ----
#' @export
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


