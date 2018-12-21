#---------------------------------------------------------------------
# Script Name:   GenotypeTableFunctions.R
# Description:   Support working with TASSEL GenotypeTables
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-21 at 15:16:56
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for
#    TASSEL classes
#--------------------------------------------------------------------

#' @title Wrapper function of TasselGenotypePhenotype class for genotype
#'    data
#' 
#' @description This function is a wrapper for the 
#'    \code{TasselGenotypePhenotype} class. It is used for storing genotype
#'    information into a class object.
#' 
#' @name readGenotypeTable
#' @rdname readGenotypeTable
#' 
#' @param path a genotype data path (e.g. \code{*.VCF, *.hmp}, etc.)
#' 
#' @importFrom rJava .jcall
#' @export
readGenotypeTable <- function(path) {
  .tasselObjectConstructor(rJava::.jcall(
    "net/maizegenetics/dna/snp/ImportUtils",
    "Lnet/maizegenetics/dna/snp/GenotypeTable;",
    "read",
    path
  )
  )
}


## Get a GenotypeTable - not exported
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

## Create Summarized Experiment from a TASSEL Genotype Table
summarizeExperimentFromGenotypeTable <- function(genotypeTable) {
  jGT <- .getTASSELClass(genotypeTable, "GenotypeTable")
  sampleDF <- sampleDataFrame(jGT)
  genomicRangesDF <- genomicRanges(jGT)
  
  genoCallIntArray <- rJava::.jcall(
    "net/maizegenetics/plugindef/GenerateRCode",
    "[I",
    "genotypeTableToDosageIntArray",
    jGT,
    FALSE
  )
  
  SummarizedExperiment(assays=matrix(genoCallIntArray,length(genomicRangesDF), byrow=TRUE), rowRanges=genomicRangesDF, colData=sampleDF)
}

## Create GWASpoly geno dataframe from SimplifiedExperiment object
GWASpolyGenoFromSummarizedExperiment <- function(SummarizedExperimentObject){
  geno <- data.frame(markerName = paste("dummy", 1:length(ranges(SummarizedExperimentObject@rowRanges)), sep = "-"), # dummy name as current summarizeExperimentFromGenotypeTable doesn't keep
                     chr = seqnames(SummarizedExperimentObject@rowRanges),
                     pos = start(ranges(SummarizedExperimentObject@rowRanges)),
                     as.data.frame(SummarizedExperimentObject@assays$data@listData) # same as assay(SummarizedExperimentObject)
  )
  colnames(geno)[4:ncol(geno)] <- as.character(SummarizedExperimentObject$Sample)
  geno
}

## Create Summarized Experiment from a TASSEL Genotype Table
snpMatrixFromGenotypeTable <- function(genotypeTable) {
  jGT <- .getTASSELClass(genotypeTable, "GenotypeTable")
  
  sampleDF <- sampleDataFrame(jGT)
  genomicRangesDF <- genomicRanges(jGT)
  
  genoCallIntArray <- rJava::.jcall(
    "net/maizegenetics/plugindef/GenerateRCode",
    "[I",
    "genotypeTableToDosageIntArray",
    jGT,
    FALSE
  )
  aMatrix <- matrix(genoCallIntArray,length(genomicRangesDF))
  colnames(aMatrix) <- sampleDF[["Sample"]]
  rownames(aMatrix) <- paste0(seqnames(genomicRangesDF),":",ranges(genomicRangesDF))
  as(aMatrix, "SnpMatrix")
}

## Create GWASpoly geno dataframe from SimplifiedExperiment object
GWASpolyGenoFromSummarizedExperiment <- function(SummarizedExperimentObject){
  geno <- data.frame(markerName = paste("dummy", 1:length(ranges(SummarizedExperimentObject@rowRanges)), sep = "-"), # dummy name as current summarizeExperimentFromGenotypeTable doesn't keep
                     chr = seqnames(SummarizedExperimentObject@rowRanges),
                     pos = start(ranges(SummarizedExperimentObject@rowRanges)),
                     as.data.frame(SummarizedExperimentObject@assays$data@listData) # same as assay(SummarizedExperimentObject)
  )
  colnames(geno)[4:ncol(geno)] <- as.character(SummarizedExperimentObject$Sample)
  geno
}