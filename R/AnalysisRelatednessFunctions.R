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


#' @title read TASSEL distance matrix object from file
#'
#' @description This function will read a TASSEL distance matrix from
#'    file and convert it into a \code{TasselDistanceMatrix} object.
#'
#' @name readDistanceMatrix
#' @rdname readDistanceMatrix
#'
#' @param file A file path of type \code{character}
#'
#' @return Returns a \code{TasselDistanceMatrix} object.
#'
#' @importFrom methods new
#' @importFrom rJava J
#'
#' @export
readDistanceMatrix <- function(file) {
    rJC <- rJava::J("net/maizegenetics/taxa/distance/ReadDistanceMatrix")
    distMatrix <- rJC$readDistanceMatrix(file)

    jTl <- distMatrix$getTaxaList()

    tl <- sapply(1:distMatrix$numberOfTaxa(), function(i) {
        jTl$taxaName(as.integer(i - 1))
    })

    methods::new(
        Class = "TasselDistanceMatrix",
        taxa = tl,
        numTaxa = distMatrix$numberOfTaxa(),
        summaryMatrix = summaryDistance(distMatrix),
        jDistMatrix = distMatrix
    )
}


#' @title Coerce matrix to TasselDistanceMatrix object
#'
#' @description Coerces an object of \code{matrix} class into an rTASSEL
#'    object of \code{TasselDistanceMatrix} class.
#'
#' @param m An object of \code{matrix} class of \eqn{m \times m} structure
#'    (e.g. a pairwise matrix). Additionally, row and column names must be
#'    the same.
#'
#' @return Returns a \code{TasselDistanceMatrix} object.
#'
#' @export
asTasselDistanceMatrix <- function(m) {

    if (nrow(m) != ncol(m)) {
        stop("Matrix object must have equal rows and columns", call. = FALSE)
    }
    if (is.null(colnames(m)) || is.null(rownames(m))) {
        stop("Matrix object must have column and row names", call. = FALSE)
    }
    if (!all(colnames(m) == rownames(m))) {
        stop("Matrix object must have the same row and column name structure", call. = FALSE)
    }

    taxa <- colnames(m)

    rJC <- rJava::J("net/maizegenetics/taxa/distance/DistanceMatrixBuilder")
    distBuilder <- rJC$getInstance(as.integer(length(taxa)))

    for(i in 1:nrow(m)) {
        for(j in 1:ncol(m)) {
            distBuilder$set(
                as.integer(j - 1),
                as.integer(i - 1),
                as.double(m[i, j])
            )
        }
        distBuilder$addTaxon(
            rJava::.jnew(
                "net/maizegenetics/taxa/Taxon",
                taxa[i]
            )
        )
    }

    distMatrix <- distBuilder$build()

    methods::new(
        Class = "TasselDistanceMatrix",
        taxa = taxa,
        numTaxa = distMatrix$numberOfTaxa(),
        summaryMatrix = summaryDistance(distMatrix),
        jDistMatrix = distMatrix
    )
}


