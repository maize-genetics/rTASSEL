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


# /// Methods (general) /////////////////////////////////////////////

## ----
#' @rdname javaRefObj
#' @aliases javaRefObj,TasselDistanceMatrix-method
#' @export
setMethod(
    f = "javaRefObj",
    signature = signature(object = "TasselDistanceMatrix"),
    definition = function(object) {
        return(object@jDistMatrix)
    }
)


## ----
#' @rdname taxaList
#' @aliases taxaList,TasselDistanceMatrix-method
#' @export
setMethod("taxaList", "TasselDistanceMatrix", function(tasObj) {
    tasObj@taxa
})



# /// Methods (coercion) ////////////////////////////////////////////

## ----
#' @title Coerce matrix from TasselDistanceMatrix class
#'
#' @description Coerces an object of class \code{TasselDistanceMatrix} to
#'    a \code{matrix} object.
#'
#' @param x An object of \code{TasselDistanceMatrix} class.
#' @param ... Additional arguments to be passed to or from methods.
#'
#' @return A \code{numeric} matrix of taxa by taxa.
#'
#' @export
as.matrix.TasselDistanceMatrix <- function(x, ...) {
    .distanceToMatrix(x@jDistMatrix, taxa = x@taxa)
}


## ----
#' @title Coerce a TasselDistanceMatrix to a dist object
#'
#' @description
#' Coerces an object of class \code{TasselDistanceMatrix} to a
#' \code{\link[stats]{dist}} object, which is the form
#' \code{\link[stats]{hclust}}, \code{\link[stats]{cmdscale}}, and most
#' other clustering functions expect.
#'
#' @details
#' A \code{dist} object holds only the lower triangle, so the diagonal of
#' a kinship matrix is dropped. Distances between a taxon and itself are
#' therefore not recoverable from the result.
#'
#' @param m An object of \code{TasselDistanceMatrix} class.
#' @param diag Should the diagonal be printed by \code{print.dist()}?
#' @param upper Should the upper triangle be printed by
#'   \code{print.dist()}?
#'
#' @return A \code{dist} object.
#'
#' @importFrom stats as.dist
#'
#' @export
as.dist.TasselDistanceMatrix <- function(m, diag = FALSE, upper = FALSE) {
    stats::as.dist(as.matrix(m), diag = diag, upper = upper)
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


