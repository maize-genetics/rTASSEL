#---------------------------------------------------------------------
# Script Name:   GenotypePhenotypeFunctions.R
# Description:   Support for TASSEL GenotypePhenotype Class
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2019-07-25 at 15:10:54
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions
#    necessary for reading in and combining phenotype and genotype
#    datasets into R and extracting data from TASSEL
#    genotype/phenotype objects.
#--------------------------------------------------------------------

#' @title Wrapper function of TasselGenotypePhenotype class for
#'    GenotypePhenotype combined data

#' @description Creates a Java GenotypePhenotype object which is used for
#'    \code{TasselGenotypePhenotype} object construction. The Java
#'    GenotypePhenotype object is created via an \code{intersect} method
#'    from TASSEL.
#'
#' @name readGenotypePhenotype
#' @rdname readGenotypePhenotype
#'
#' @param genoPathOrObj a path to a genotype file (e.g. VCF, hmp, etc.) or
#'    TASSEL Genotype Obj
#' @param phenoPathDFOrObj a path, a data frame of phenotypic data, or TASSEL
#'    Phenotype Obj
#' @param ... Additional parameters to be sent to the function. Currently,
#'    if an R data frame object is passed, additional parameters will also
#'    need to be entered for this process. These parameters are derived from
#'    the \code{readPhenotypeFromDataFrame} function. Mainly, \code{taxaID}
#'    is required. If you would like to specify depth retention and position
#'    sorting in Genotype Tables from a path, indicate them here. See
#'    \code{readGenotypeTableFromPath()} for more detail.
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @importFrom rJava J
#' @importFrom rJava is.jnull
#' @importFrom rJava .jinstanceof
#' @export
readGenotypePhenotype <- function(genoPathOrObj, phenoPathDFOrObj, ...) {
    genoObj <- getGenotypeTable(genoPathOrObj)
    if(rJava::is.jnull(genoObj)) {
        genoObj <- getGenotypeTable(
            readGenotypeTableFromPath(genoPathOrObj, keepDepth = FALSE, sortPositions = FALSE)
        )
    }
    phenoObj <- getPhenotypeTable(phenoPathDFOrObj)
    if(rJava::is.jnull(phenoObj) & is.data.frame(phenoPathDFOrObj)) {
        phenoObj <- getPhenotypeTable(
            readPhenotypeFromDataFrame(phenoPathDFOrObj, ...)
        )
    } else if(rJava::is.jnull(phenoObj) & !is.data.frame(phenoPathDFOrObj)) {
        phenoObj <- rJava::new(
            rJava::J("net/maizegenetics/phenotype/PhenotypeBuilder")
        )$fromFile(phenoPathDFOrObj)$build()$get(0L)
    }

    t <- .tasselObjectConstructor(
        combineTasselGenotypePhenotype(
            genotypeTable = genoObj,
            phenotype = phenoObj
        )
    )
    return(t)
}


## Combine TASSEL GenotypeTable and Phenotype to GenotypePhenotype class - not exported (house keeping)
combineTasselGenotypePhenotype <- function(genotypeTable, phenotype) {
    genotypeTable <- getGenotypeTable(genotypeTable)
    phenotype <- getPhenotypeTable(phenotype)
    new(J("net.maizegenetics.phenotype.GenotypePhenotypeBuilder"))$
        genotype(genotypeTable)$phenotype(phenotype)$intersect()$build()
}


## Get GenotypePhenotype - not exported - not exported (house keeping)
getGenotypePhenotype <- function(jtsObject) {
    if(is(jtsObject, "TasselGenotypePhenotype")) {
        jtsObject <- jtsObject@jTasselObj
    }
    if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
    if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
        return(jtsObject)
    } else {
        return(rJava::.jnull())
    }
}
