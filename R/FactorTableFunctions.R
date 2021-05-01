#---------------------------------------------------------------------
# Description:   Support working with TASSEL FactorTables
# Author:        Brandon Monier & Terry Casstevens
# Created:       2021-04-30 at 10:07:53
# Last Modified: 2021-04-30 at 10:07:53
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions
#    necessary for reading in haplotype/factor datasets into R and
#    extracting data from TASSEL 6 objects.
#--------------------------------------------------------------------

#' @title Wrapper function of FactorTable class for haplotype data data
#'
#' @description This function is a wrapper for the
#'    \code{FactorTable} class. It is used for storing genotype
#'    information into a class object.
#'
#' @name readFactorTable
#' @rdname readFactorTable
#'
#' @param file A haplotype or factor data path (e.g. \code{*.vcf}).
#' @param verbose Show messages to console? Defaults to \code{TRUE}.
#'
#' @return Returns an object of \code{FactorTable} class.
#'
#' @importFrom rJava .jnew
#' @importFrom methods new
#' @export
readFactorTable <- function(file, verbose = TRUE) {

    if (!file.exists(file)) {
        stop("Cannot open file ", file, ": No such file or directory")
    }

    message("Reading VCF file...")
    factorBuilder <- rJava::.jnew("net/maizegenetics/dna/factor/io/BuilderFromHaplotypeVCF")

    methods::new(
        Class = "FactorTable",
        jFactorTable = factorBuilder$read(file)
    )
}


#' @title Get reference ranges
#'
#' @description Returns reference range intervals from a \code{FactorTable}
#'    object.
#'
#' @name refRanges
#' @rdname refRanges
#'
#' @param factorTable An object of clase \code{FactorTable}.
#'
#' @return Returns a \code{GRanges} object of reference range data.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom rJava .jevalArray
#' @importFrom rJava .jnew
#' @importFrom S4Vectors Rle
#' @export
refRanges <- function(factorTable) {
    jRC <- rJava::.jnew("net/maizegenetics/plugindef/GenerateRCodeT6")
    rr <- jRC$getRefRanges(factorTable@jFactorTable)
    rr <- lapply(rJava::.jevalArray(rr), rJava::.jevalArray)

    GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(rr[[1]]),
        ranges = IRanges::IRanges(
            start = as.numeric(rr[[3]]),
            end = as.numeric(rr[[4]])
        )
    )
}


#' @title Get taxa
#'
#' @description Returns taxa IDs from a \code{FactorTable} object.
#'
#' @name taxa
#' @rdname taxa
#'
#' @param factorTable An object of clase \code{FactorTable}.
#'
#' @return Returns a vector of taxa of type \code{character}.
#'
#' @importFrom rJava .jnew
#' @export
taxa <- function(factorTable) {
    jRC <- rJava::.jnew("net/maizegenetics/plugindef/GenerateRCodeT6")
    jRC$getTaxaArray(factorTable@jFactorTable)
}


#' @title Get hap ID matrix
#'
#' @importFrom rJava .jnew
#'
#' @description returns a matrix of hap ID integers.
factorTableHapMatrix <- function(factorTable) {
    jRC <- rJava::.jnew("net/maizegenetics/plugindef/GenerateRCodeT6")
    intArray <- jRC$getHapArray(factorTable@jFactorTable)
    m <- lapply(intArray, rJava::.jevalArray)
    m <- simplify2array(m)
    rownames(m) <- taxa(factorTable)
    colnames(m) <- 1:factorTable@jFactorTable$numFactors()
    return(m)
}


#' @title Get haplotype ID data frame
#'
#' @description Returns a \code{data.frame} of haplotype and taxa IDs from a
#'    \code{FactorTable} object.
#'
#' @name factorTableDF
#' @rdname factorTableDF
#'
#' @param factorTable An object of clase \code{FactorTable}.
#'
#' @return Returns a \code{data.frame} of haplotype and taxa IDs.
#'
#' @export
factorTableDF <- function(factorTable) {
    taxa <- data.frame(taxa = taxa(factorTable))
    m <- factorTableHapMatrix(factorTable)
    m <- as.data.frame(m, row.names = NULL)
    m[,] <- lapply(m[,], factor)

    m <- cbind(taxa, m)
    rownames(m) <- NULL

    return(m)
}


#' @title Create Summarized Experiment from a TASSEL Genotype Table
#'
#' @description This function will generate an object of
#'    \code{SummarizedExperiment} class for marker data derived from a
#'    \code{FactorTable} class object.
#'
#' @name getSumExp
#' @rdname getSumExp
#'
#' @param factorTable An object of clase \code{FactorTable}.
#'
#' @return Returns a \code{SummarizedExperiment} of TASSEL genotype data.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
getSumExp <- function(factorTable) {
    SummarizedExperiment::SummarizedExperiment(
        assays = list(hapID = t(factorTableHapMatrix(factorTable))),
        rowRanges = refRanges(factorTable)
    )
}


