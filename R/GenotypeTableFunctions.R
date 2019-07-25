#---------------------------------------------------------------------
# Script Name:   GenotypeTableFunctions.R
# Description:   Support working with TASSEL GenotypeTables
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2019-07-25 at 14:55:30
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions
#    necessary for reading in genotype datasets into R and
#    extracting data from TASSEL genotype objects.
#--------------------------------------------------------------------

#' @title Wrapper function of TasselGenotypePhenotype class for genotype
#'    data
#'
#' @description This function is a wrapper for the
#'    \code{TasselGenotypePhenotype} class. It is used for storing genotype
#'    information into a class object.
#'
#' @name readGenotypeTableFromPath
#' @rdname readGenotypeTableFromPath
#'
#' @param path A genotype data path (e.g. \code{*.VCF, *.hmp}, etc.).
#' @param keepDepth Should depth be kept? Defaults to \code{FALSE}.
#' @param sortPositions Should positions be sorted? Defaults to \code{FALSE}.
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @importFrom rJava J
#' @importFrom rJava %instanceof%
#' @export
readGenotypeTableFromPath <- function(path, keepDepth = FALSE, sortPositions = FALSE) {
    if (!file.exists(path)) {
        stop("Cannot open file ", path, ": No such file or directory")
    }

    jrc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    return(
        .tasselObjectConstructor(
            jrc$read(path, keepDepth, sortPositions)
        )
    )
}


#' @title Create Summarized Experiment from a TASSEL Genotype Table
#'
#' @description This function will generate an object of
#'    \code{SummarizedExperiment} class for marker data derived from a
#'    \code{TasselGenotypePhenotype} class object.
#'
#' @name getSumExpFromGenotypeTable
#' @rdname getSumExpFromGenotypeTable
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#' @param useRef Use reference alleles or major alleles at sites. If
#'    \code{FALSE}, major alleles will be used.
#'
#' @return Returns a \code{SummarizedExperiment} of TASSEL genotype data.
#'
#' @importFrom rJava .jcall
#' @importFrom rJava is.jnull
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
getSumExpFromGenotypeTable <- function(tasObj, useRef = FALSE) {
    jGT <- getGenotypeTable(tasObj)

    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    if (rJava::is.jnull(jGT)) {
        stop("TASSEL genotype table not detected.")
    }

    sampleDF <- sampleDataFrame(jGT)
    genomicRangesDF <- genomicRanges(jGT)

    jc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")

    genoCallIntArray <- jc$genotypeTableToDosageIntArray(
        jGT,
        useRef
    )

    SummarizedExperiment::SummarizedExperiment(
        assays = matrix(
            genoCallIntArray,
            length(genomicRangesDF),
            byrow = TRUE
        ),
        rowRanges = genomicRangesDF,
        colData = sampleDF
    )
}


## Get a GenotypeTable - not exported (house keeping)
getGenotypeTable <- function(jtsObject) {
    if(is(jtsObject, "TasselGenotypePhenotype")) {
        return(jtsObject@jGenotypeTable)
    }
    if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
    if(jtsObject %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
        return(jtsObject)
    } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
        return(jtsObject$genotypeTable())
    } else {
        return(rJava::.jnull())
    }
}
