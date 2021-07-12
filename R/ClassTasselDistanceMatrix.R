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
        cat("A TasselDistanceMatrix Object\n")
        cat("  n x n: ", object@jDistMatrix$numberOfTaxa)
    }
)


