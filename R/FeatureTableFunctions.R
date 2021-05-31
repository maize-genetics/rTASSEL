#---------------------------------------------------------------------
# Description:   Support working with TASSEL FeatureTables
# Author:        Brandon Monier & Terry Casstevens
# Created:       2021-04-30 at 10:07:53
# Last Modified: 2021-05-28 at 12:59:55
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions
#    necessary for reading in haplotype/factor datasets into R and
#    extracting data from TASSEL 6 objects.
#--------------------------------------------------------------------

#' @title Read feature information from \code{HVCF} files
#'
#' @description This function is a wrapper for the
#'    \code{FeatureTable} class. It is used for storing genotype feature
#'    information into a class object.
#'
#' @name readFeatureTableFromVCF
#' @rdname readFeatureTableFromVCF
#'
#' @param file A haplotype or feature data path (e.g. \code{*.vcf}).
#' @param verbose Show messages to console? Defaults to \code{TRUE}.
#'
#' @return Returns an object of \code{FeatureTable} class.
#'
#' @importFrom methods new
#' @importFrom rJava .jnew
#' @export
readFeatureTableFromVCF <- function(file, verbose = TRUE) {

    if (!file.exists(file)) {
        stop("Cannot open file ", file, ": No such file or directory")
    }

    if (verbose) message("Reading HVCF file...")
    featureBuilder <- rJava::.jnew("net/maizegenetics/dna/factor/io/BuilderFromHaplotypeVCF")

    methods::new(
        Class = "FeatureTable",
        jFeatureTable = featureBuilder$read(file)
    )
}


#' @title Read feature information from a BrAPI web service
#'
#' @description Returns PHG feature information hosted on a BrAPI server.
#'
#' @name readFeatureTableFromBrapi
#' @rdname readFeatureTableFromBrapi
#'
#' @param brapiObj An BrAPI connection object of \code{BrapiCon} class.
#' @param method A BrAPI method ID, currently linked to PHG methods.
#' @param verbose Show messages to console? Defaults to \code{TRUE}.
#'
#' @return Returns an object of \code{FeatureTable} class.
#'
#' @importFrom methods new
#' @importFrom rJava .jnew
#' @importFrom rPHG brapiURL
#' @export
readFeatureTableFromBrapi <- function(brapiObj, method, verbose = TRUE) {

    if (class(brapiObj) != "BrapiCon") {
        stop("`tasObj` must be of class `BrapiCon`", call. = FALSE)
    }

    rJC <- rJava::.jnew(
        "net/maizegenetics/analysis/brapi/BrAPIConnection",
        rPHG::brapiURL(brapiObj)
    )

    if (verbose) message("Downloading BrAPI feature information...")
    methods::new(
        Class = "FeatureTable",
        jFeatureTable = rJC$getCalls(method)
    )
}


#' @title Get reference ranges
#'
#' @description Returns reference range intervals from a \code{FeatureTable}
#'    object.
#'
#' @name refRanges
#' @rdname refRanges
#'
#' @param featureTable An object of clase \code{FeatureTable}.
#'
#' @return Returns a \code{GRanges} object of reference range data.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom rJava .jevalArray
#' @importFrom rJava .jnew
#' @importFrom S4Vectors Rle
#' @export
refRanges <- function(featureTable) {
    jRC <- rJava::.jnew("net/maizegenetics/plugindef/RCodeHelpers")
    rr <- jRC$getRefRanges(featureTable@jFeatureTable)
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
#' @description Returns taxa IDs from a \code{FeatureTable} object.
#'
#' @name taxa
#' @rdname taxa
#'
#' @param featureTable An object of clase \code{FeatureTable}.
#'
#' @return Returns a vector of taxa of type \code{character}.
#'
#' @importFrom rJava .jnew
#' @export
taxa <- function(featureTable) {
    jRC <- rJava::.jnew("net/maizegenetics/plugindef/RCodeHelpers")
    jRC$getTaxaArray(featureTable@jFeatureTable)
}


#' @title Get hap ID matrix
#'
#' @importFrom rJava .jevalArray
#' @importFrom rJava .jnew
#'
#' @description returns a matrix of hap ID integers.
featureTableHapMatrix <- function(featureTable) {
    jRC <- rJava::.jnew("net/maizegenetics/plugindef/RCodeHelpers")
    intArray <- jRC$getHapArray(featureTable@jFeatureTable)
    m <- lapply(intArray, rJava::.jevalArray)
    m <- simplify2array(m)
    rownames(m) <- taxa(featureTable)
    colnames(m) <- 1:featureTable@jFeatureTable$numFeatures()
    return(m)
}


#' @title Get haplotype ID data frame
#'
#' @description Returns a \code{data.frame} of haplotype and taxa IDs from a
#'    \code{FeatureTable} object.
#'
#' @name featureTableDF
#' @rdname featureTableDF
#'
#' @param featureTable An object of clase \code{FeatureTable}.
#'
#' @return Returns a \code{data.frame} of haplotype and taxa IDs.
#'
#' @export
featureTableDF <- function(featureTable) {
    taxa <- data.frame(taxa = taxa(featureTable))
    m <- featureTableHapMatrix(featureTable)
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
#'    \code{FeatureTable} class object.
#'
#' @name getSumExp
#' @rdname getSumExp
#'
#' @param featureTable An object of clase \code{FeatureTable}.
#'
#' @return Returns a \code{SummarizedExperiment} of TASSEL genotype data.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
getSumExp <- function(featureTable) {
    SummarizedExperiment::SummarizedExperiment(
        assays = list(hapID = t(featureTableHapMatrix(featureTable))),
        rowRanges = refRanges(featureTable)
    )
}


