## ----
#' @title TasselDistanceMatrix Class
#'
#' @description Class \code{TasselDistanceMatrix} defines a \code{rTASSEL}
#'    Class for storing TASSEL genotype and phenotype objects.
#'
#' @name TasselDistanceMatrix-class
#' @rdname TasselDistanceMatrix-class
#' @exportClass TasselDistanceMatrix
#' @importFrom BiocGenerics colnames rownames ncol nrow
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


#' @title Show method TasselDistanceMatrix objects
#'
#' @description Prints out information related taxa, positions, genotype, and
#'    phenotype information.
#'
#' @param object a \code{TasselGenotypePhenotype} class object
#'
#' @rdname TasselDistanceMatrix-class
#' @aliases show,TasselDistanceMatrix-method
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
#' @param x An object of \code{TasselDistanceMatrix} class.
#' @param ... Additional arguments to be passed to or from methods.
#'
#' @export
as.matrix.TasselDistanceMatrix <- function(x, ...) {
    tmp1 <- unlist(strsplit(x@jDistMatrix$toStringTabDelim(), split = "\n"))
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
#'
#' @param x An object of class \code{TasselDistanceMatrix}.
setMethod(
    f = "dim",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        c(x@numTaxa, x@numTaxa)
    }
)

#' @title Column names
#'
#' @description Get column names of a \code{TasselDistanceMatrix} object.
#'
#' @param x An object of class \code{TasselDistanceMatrix}.
#'
#' @export
setMethod(
    f = "colnames",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        x@taxa
    }
)


#' @title Row names
#'
#' @description Get row names of a \code{TasselDistanceMatrix} object.
#'
#' @param x An object of class \code{TasselDistanceMatrix}.
#'
#' @export
setMethod(
    f = "rownames",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        x@taxa
    }
)


#' @title Number of columns
#'
#' @description Get number of columns of a \code{TasselDistanceMatrix} object.
#'
#' @param x An object of class \code{TasselDistanceMatrix}.
#'
#' @export
setMethod(
    f = "ncol",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        x@numTaxa
    }
)


#' @title Number of rows
#'
#' @description Get number of rows of a \code{TasselDistanceMatrix} object.
#'
#' @param x An object of class \code{TasselDistanceMatrix}.
#'
#' @export
setMethod(
    f = "nrow",
    signature = "TasselDistanceMatrix",
    definition = function(x) {
        x@numTaxa
    }
)


