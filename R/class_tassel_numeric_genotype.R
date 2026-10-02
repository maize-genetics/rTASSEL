## ----
#' @title
#' TasselNumericGenotype Class Definition
#'
#' @description
#' Defines the `TasselNumericGenotype` class, which extends the
#' `TasselGenotype` class. This class is used to represent numeric
#' genotype data in the TASSEL 5 framework.
#'
#' @include class_tassel_genotype.R
#'
#' @name TasselNumericGenotype-class
#' @rdname TasselNumericGenotype-class
#' @exportClass TasselNumericGenotype
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
#' @param object
#' An object of class \code{TasselNumericGenotype}.
#'
#' @docType methods
#' @rdname TasselNumericGenotype-class
#' @aliases show,TasselNumericGenotype-method
setMethod("show", "TasselNumericGenotype", function(object) {
    fgs <- formatNumGtStrings(object@jRefObj, nTaxa = 5, nSites = 5)
    printGtDisp(
        fgs       = fgs,
        nTaxa     = object@jRefObj$numberOfTaxa(),
        nSites    = object@jRefObj$numberOfSites(),
        jMem      = object@jMemAddress,
        className = "TasselNumericGenotype"
    )
})



# /// Methods (coercion) ////////////////////////////////////////////

## ----
#' @title Coerce numeric genotype data to an R matrix
#'
#' @description
#' Converts the numeric genotype table held by a
#' \code{TasselNumericGenotype} object into a matrix of reference
#' probabilities, with taxa as rows and sites as columns.
#'
#' @details
#' A numeric genotype table holds a probability rather than a discrete
#' call, so there is no dosage to report and the \code{type} argument of
#' \code{\link{as.matrix.TasselGenotype}} does not apply.
#'
#' TASSEL 5 reads reference probabilities one cell at a time, so this
#' costs a Java call per cell and is far slower than the dosage matrix of
#' a comparably sized allele-based table. Filter the table down to the
#' taxa and sites of interest before materializing it.
#'
#' @param x A \code{TasselNumericGenotype} object.
#' @param ... Additional arguments to be passed to or from methods.
#'
#' @return A \code{numeric} matrix of taxa (rows) by sites (columns).
#'
#' @examples
#' \dontrun{
#' numGtPath <- system.file("extdata", "numeric_genotype.txt", package = "rTASSEL")
#'
#' readGenotype(numGtPath) |> as.matrix()
#' }
#'
#' @export
as.matrix.TasselNumericGenotype <- function(x, ...) {
    if (!x@jRefObj$hasReferenceProbablity()) {
        rlang::abort(c(
            "`x` does not contain reference probabilities",
            "i" = "Only numeric genotype tables can be coerced this way"
        ))
    }

    .refProbMatrix(
        jGt       = x@jRefObj,
        taxa      = taxaList(x),
        siteNames = positionList(x)$Name
    )
}


