## ----
#' @name TasselGenotype-class
#' @rdname TasselGenotype-class
#' @exportClass TasselGenotype
setClass(
    Class = "TasselGenotype",
    slots = c(
        dispData    = "list",
        jRefObj     = "jobjRef",
        jMemAddress = "character",
        jClass      = "character"
    )
)


## ----
# targets:
#   * path
#   * matrix
#   * data frame
#' @title Helper function to build a TasselGenotype object
#' @export
readGenotype <- function(x, sortPositions = FALSE, keepDepth = FALSE) {
    if (is.character(x)) {
        xNorm <- normalizePath(x, mustWork = FALSE)
        if (!file.exists(xNorm)) {
            rlang::abort("The input path is not a valid file")
        }

        readGenotypeFromPath(x, sortPositions, keepDepth)
    } else if (is.matrix(x)) {
        readGenotypeFromRMatrix(x)
    } else {
        rlang::abort("Unsupported data type")
    }
}



# /// Methods (show) ////////////////////////////////////////////////

## ----
#' @export
setMethod("show", "TasselGenotype", function(object) {
    printGtDisp(
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
    signature = signature(object = "TasselGenotype"),
    definition = function(object) {
        return(object@jRefObj)
    }
)



