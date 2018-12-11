#---------------------------------------------------------------------
# Script Name:   PullFunctions.R
# Description:   Various tests with rJava
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-03 at 17:58:46
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for 
#    TASSEL classes
#--------------------------------------------------------------------

source("R/AllClasses.R")


## Constructor for GenotypeTable class object
sampleDataFrame <- function(jtsGenoTableOrTaxaList) {
  if(is(jtsGenoTableOrTaxaList,"GenotypeTable")) {
    jtsTL <- taxa(jtsGenoTableOrTaxaList)@jtsTaxaList
  } else if(is(jtsGenoTableOrTaxaList,"TaxaList")) {
    jtsTL <- jtsGenoTableOrTaxaList@jtsTaxaList
  } else if(jtsGenoTableOrTaxaList %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
    jtsTL <- jtsGenoTableOrTaxaList$taxa()
  } else if(jtsGenoTableOrTaxaList %instanceof% "net.maizegenetics.taxa.TaxaList") {
    jtsTL <- jtsGenoTableOrTaxaList
  } else {
    stop("Object is not of \"TaxaList\" class")
  }
  
  #This could be faster
  # taxaArray <- rJava::.jcall(
  #   "net/maizegenetics/plugindef/GenerateRCode",
  #   "[S",
  #   "genotypeTableToSampleNameArray",
  #   test@jtsGenotypeTable
  # )
  taxaArray <- c()
  for(i in 1:jtsTL$size()) {
    #why do I have to do -1L
    taxaArray[i] = toString(jtsTL$taxaName(i-1L))
  }
  colData <- data.frame(Sample=taxaArray, row.names = taxaArray)
  colData
}


## Constructor for GRanges (GenomicRanges) class object
genomicRanges <- function(genoTable) {
    if(is(genoTable,"GenotypeTable")) {
        jtsPL <- positions(genoTable)@jtsPositionList
    } else if(genoTable %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
      jtsPL <- genoTable$positions()
    } else {
        stop("Object is not of \"GenotypeTable\" class")
    }
    
  # genoPositionVector <- rJava::.jcall(
  #   "net/maizegenetics/plugindef/GenerateRCode",
  #   "Lnet/maizegenetics/plugindef/GenerateRCode$PositionVectors;",
  #   "genotypeTableToPositionListOfArrays",
  #   test@jtsGenotypeTable
  # )
  
    numSite <- as.numeric(jtsPL$numberOfSites())
    physPos <- jtsPL$physicalPositions()
    
    cat("Extracting chromosome names for each postion...\n")
    cat("...is there a quicker way to get this? (~ Brandon)\n")
    chrName <- lapply(seq_len(numSite), function(pos) {
        jtsPL$chromosomeName(as.integer(pos - 1))
    })
    chrName <- unlist(chrName)
    
    gr2 <- GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(chrName),
        ranges = IRanges::IRanges(start = physPos, end = physPos)
    )
    return(gr2)
}

## Create Summarized Experiment from a TASSEL Genotype Table
summarizeExperimentFromGenotypeTable <- function(genotypeTable) {
  if(is(genotypeTable,"GenotypeTable")) {
    jGT <- genotypeTable@jtsGenotypeTable
  } else if(genotypeTable %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
    jGT <- genotypeTable
  } else {
    stop("Object is not of \"GenotypeTable\" class")
  }
  
  sampleDF <- sampleDataFrame(jGT)
  genomicRangesDF <- genomicRanges(jGT)
  
  genoCallIntArray <- rJava::.jcall(
    "net/maizegenetics/plugindef/GenerateRCode",
    "[I",
    "genotypeTableToDosageIntArray",
    jGT
  )
  
 SummarizedExperiment(assays=matrix(genoCallIntArray,length(genomicRangesDF)), rowRanges=genomicRangesDF, colData=sampleDF)
}


