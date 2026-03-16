## ----
#' @title TasselGenotype Class
#' @description An S4 class to represent a Tassel Genotype object.
#'
#' @slot jRefObj
#' A reference to a Java object (`jobjRef`) associated with the genotype.
#' @slot jMemAddress
#' A character string representing the memory address of the Java object.
#' @slot jClass
#' A character string representing the Java class of the object.
#'
#' @details
#' This class is designed to interface with TASSEL 5 for genotype
#' data management and analysis. It provides a structure to store and
#' interact with Java objects used in TASSEL 5.
#'
#' @name TasselGenotype-class
#' @rdname TasselGenotype-class
#' @exportClass TasselGenotype
setClass(
    Class = "TasselGenotype",
    slots = c(
        jRefObj     = "jobjRef",
        jMemAddress = "character",
        jClass      = "character"
    )
)


## ----
#' @title
#' Read Genotype Data
#'
#' @description
#' This function reads genotype data from a file path or an R matrix.
#' It supports optional sorting of positions and retaining depth
#' information.
#'
#' @details
#' \itemize{
#'   \item If \code{x} is a character string:
#'     \itemize{
#'       \item The function checks if the file exists.
#'       \item Reads the genotype data from the file path using
#'       \code{readGenotypeFromPath}.
#'     }
#'   \item If \code{x} is a matrix:
#'     \itemize{
#'       \item The function processes the genotype data using
#'       \code{readGenotypeFromRMatrix}.
#'     }
#'   \item If \code{x} is neither a character string nor a matrix:
#'     \itemize{
#'       \item An error is raised.
#'     }
#' }
#'
#' @param x
#' A character string representing the file path to the genotype data
#' or a matrix containing genotype data.
#' @param sortPositions
#' A logical value indicating whether to sort positions in the
#' genotype data. Default is \code{FALSE}.
#' @param keepDepth
#' A logical value indicating whether to retain depth information in
#' the genotype data. Default is \code{FALSE}.
#'
#' @examples
#' \dontrun{
#' # Read genotype data from a file
#' readGenotype("path/to/genotype/file.txt", sortPositions = TRUE, keepDepth = TRUE)
#'
#' # Read genotype data from a matrix
#' genotypeMatrix <- matrix(data = ..., nrow = ..., ncol = ...)
#' readGenotype(genotypeMatrix)
#' }
#'
#' @return
#' A processed genotype object based on the input data.
#'
#' @export
readGenotype <- function(x, sortPositions = FALSE, keepDepth = FALSE) {
    if (is.character(x)) {
        xNorm <- normalizePath(x, mustWork = FALSE)
        if (!file.exists(xNorm)) {
            rlang::abort("The input path is not a valid file")
        }

        readGenotypeFromPath(xNorm, sortPositions, keepDepth)
    } else if (is.matrix(x)) {
        readNumericGenotypeFromRMatrix(x, asTGP = FALSE)
    } else {
        rlang::abort("Unsupported data type")
    }
}



# /// Methods (show) ////////////////////////////////////////////////

## ----
#' @title
#' Display TasselGenotype Object
#'
#' @description
#' This method is used to display a summary of a `TasselGenotype`
#' object. It prints genotype display information, including the
#' number of taxa, number of sites, and memory address of the Java
#' object.
#'
#' @param object
#' An object of class `TasselGenotype`.
#'
#' @method show TasselGenotype
setMethod("show", "TasselGenotype", function(object) {
    fgs <- formatGtStrings(object@jRefObj)
    printGtDisp(
        fgs    = fgs,
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



# /// Methods (summary) /////////////////////////////////////////////

## ----
#' @rdname siteSummary
#' @aliases siteSummary,TasselGenotype-method
#' @export
setMethod("siteSummary", "TasselGenotype", function(tasObj) {
    .runGenotypeSummary(tasObj@jRefObj, doSite = TRUE)
})


## ----
#' @rdname taxaSummary
#' @aliases taxaSummary,TasselGenotype-method
#' @export
setMethod("taxaSummary", "TasselGenotype", function(tasObj) {
    .runGenotypeSummary(tasObj@jRefObj, doTaxa = TRUE)
})



# /// Bracket Method /////////////////////////////////////////////////

## ----
#' @title Subset a TasselGenotype
#'
#' @description
#' Matrix-style subsetting for \code{TasselGenotype} objects using
#' the \code{gt[taxa, sites]} syntax.
#'
#' @param x A \code{TasselGenotype} object.
#' @param i Taxa selector: a character vector of IDs, a
#'   \code{\linkS4class{TaxaSelector}}, or missing.
#' @param j Site selector: an integer vector of 0-based indices, a
#'   character vector of site names, a
#'   \code{\linkS4class{SiteSelector}}, or missing.
#' @param ... Ignored.
#' @param drop Ignored.
#'
#' @return A new \code{TasselGenotype} (or subclass) containing the
#'   selected taxa and/or sites.
#'
#' @examples
#' \dontrun{
#' gt[taxa("B73", "Mo17"), ]
#' gt[, sites(0:999)]
#' gt[, sitesWhere(maf >= 0.05)]
#' gt[taxa("B73"), region("chr1", 1e6, 2e6)]
#' gt[, !sites(0:9)]
#' }
#'
#' @rdname TasselGenotype-class
#' @aliases [,TasselGenotype,ANY,ANY-method
setMethod("[", "TasselGenotype", function(x, i, j, ..., drop = FALSE) {
    jGt <- x@jRefObj
    if (!missing(i)) jGt <- applyTaxaSelector(jGt, i)
    if (!missing(j)) jGt <- applySiteSelector(jGt, j)
    newTasselGenotype(jGt, x)
})
