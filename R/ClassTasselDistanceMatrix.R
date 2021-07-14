#--------------------------------------------------------------------
# Script Name:   ClassTasselDistanceMatrix.R
# Description:   TasselDistanceMatrix class and methods
# Author:        Brandon Monier
# Created:       2021-07-12 at 15:03:58
# Last Modified: 2021-07-12 at 15:05:38
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# TasselDistanceMatrix class, constructors, and methods
#--------------------------------------------------------------------

#' @title TasselDistanceMatrix Class
#'
#' @description Class \code{TasselDistanceMatrix} defines a \code{rTASSEL}
#'    Class for storing TASSEL genotype and phenotype objects.
#'
#' @name TasselDistanceMatrix-class
#' @rdname TasselDistanceMatrix-class
#' @exportClass TasselDistanceMatrix
setClass(
    Class = "TasselDistanceMatrix",
    representation = representation(
        name = "character",
        taxa = "character",
        numTaxa = "numeric",
        summaryMatrix = "matrix",
        jDistMatrix = "jobjRef"
    )
)


#' @title Show method TasselGenotypePhenotype objects
#'
#' @description Prints out information related taxa, positions, genotype, and
#'    phenotype information.
#'
#' @param object a \code{TasselGenotypePhenotype} class object
#'
#' @rdname TasselGenotypePhenotype-class
#' @aliases show,TasselGenotypePhenotype-method
#'
#' @importFrom rJava .jnull
setMethod(
    f = "show",
    signature = "TasselDistanceMatrix",
    definition = function(object) {
        m <- object@jDistMatrix$numberOfTaxa()
        s <- object@summaryMatrix
        cat("A TasselDistanceMatrix Object of", m, "x", m, "elements:")
        cat("\n\n")
        for (i in seq_len(nrow(s))) {
            cat(" ", object@summaryMatrix[i, ])
            cat("\n")
        }
    }
)


# @TODO - Optimize this code to use 2d array methods
#' @title Coerce matrix from TasselDistanceMatrix class
#'
#' @description Coerces an object of class \code{TasselDistanceMatrix} to
#'    a \code{matrix} object.
#'
#' @param object An object of \code{TasselDistanceMatrix} class.
#'
#' @export
as.matrix.TasselDistanceMatrix <- function(object) {
    tmp1 <- unlist(strsplit(object@jDistMatrix$toStringTabDelim(), split = "\n"))
    tmp2 <- strsplit(tmp1, split = "\t")
    tmp3 <- t(simplify2array(tmp2))
    colnames(tmp3) <- as.character(unlist(tmp3[1, ]))
    tmp3 <- tmp3[-1, ]
    matRow <- tmp3[, 1]
    tmp3 <- tmp3[, -1]
    tmp3 <- apply(tmp3, 2, as.numeric)
    rownames(tmp3) <- matRow
    return(tmp3)
}


#' @title Get dimensions of TasselDistanceMatrix object
#'
#' @description Retrieves dimensions of a \code{TasselDistanceMatrix} object
#'    as a vector (e.g. \code{c(10, 10)}).
setMethod(
    f = "dim",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        c(x@numTaxa, x@numTaxa)
    }
)

#' @export
setMethod(
    f = "colnames",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        x@taxa
    }
)


#' @export
setMethod(
    f = "rownames",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        x@taxa
    }
)


#' @export
setMethod(
    f = "ncol",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        x@numTaxa
    }
)


#' @export
setMethod(
    f = "nrow",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        x@numTaxa
    }
)







