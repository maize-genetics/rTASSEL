#--------------------------------------------------------------------
# Script Name:   AnalysisRelatednessFunctions.R
# Description:   Functions for TASSEL relatedness analyses
# Author:        Brandon Monier
# Created:       2019-04-04 at 21:31:09
# Last Modified: 2019-04-04 at 22:25:19
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions necessary
#    for relatedness-based methods in TASSEL.
#--------------------------------------------------------------------

#' @title Create a TASSEL kinship matrix
#'
#' @description This function will calculate a kinship matrix using
#'    TASSEL's \code{Kinship} plugin and respective parameters.
#'
#' @name kinshipMatrix
#' @rdname kinshipMatrix
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#' @param method A Kinship method.
#' @param maxAlleles Maximum number of alleles.
#' @param algorithmVariation Algorithm variation.
#'
#' @return Returns a Java pointer of a TASSEL kinship matrix object.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#' @export
kinshipMatrix <- function(tasObj,
                          method = "Centered_IBS",
                          maxAlleles = 6,
                          algorithmVariation = "Observed_Allele_Freq") {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    # Create kinship plugin
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.distance.KinshipPlugin"),
        rJava::.jnull(),
        FALSE
    )
    plugin$setParameter("method", toString(method))
    plugin$setParameter("maxAlleles", toString(maxAlleles))
    plugin$setParameter("algorithmVariation", toString(algorithmVariation))
    plugin$runPlugin(jGenoTable)
}


#' @title Convert TASSEL kinship matrix object to an R matrix class
#'
#' @description This function will take a TASSEL kinship object and convert
#'    it into an R \code{matrix} object.
#'
#' @name kinshipToRMatrix
#' @rdname kinshipToRMatrix
#'
#' @param kinJobj A TASSEL kinship object.
#'
#' @return Returns an R \code{matrix} object.
#'
#' @export
kinshipToRMatrix <- function(kinJobj) {
    tmp1 <- unlist(strsplit(kinJobj$toStringTabDelim(), split = "\n"))
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


#' @title Create a TASSEL distance matrix
#'
#' @description This function will calculate a distance matrix using
#'    TASSEL's \code{DistanceMatrix} plugin.
#'
#' @name distanceMatrix
#' @rdname distanceMatrix
#'
#' @param tasObj an rTASSEL \code{TasselGenotypePhenotype} object
#'
#' @return Returns a Java pointer of a TASSEL distance matrix object.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#' @export
distanceMatrix <- function(tasObj) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.distance.DistanceMatrixPlugin"),
        rJava::.jnull(),
        FALSE
    )
    plugin$getDistanceMatrix(jGenoTable)
}


#' @title Convert TASSEL distance matrix object to an R matrix class
#'
#' @description This function will take a TASSEL distance matrix object and
#'    convert it into an R \code{matrix} object.
#'
#' @name distanceToRMatrix
#' @rdname distanceToRMatrix
#'
#' @param distJobj A TASSEL distance matrix object.
#'
#' @return Returns an R \code{matrix} object.
#'
#' @export
distanceToRMatrix <- function(distJobj) {
    tmp1 <- unlist(strsplit(distJobj$toStringTabDelim(), split = "\n"))
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
