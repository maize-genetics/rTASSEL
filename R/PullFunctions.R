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

  taxaArray <- J("net/maizegenetics/plugindef/GenerateRCode")$genotypeTableToSampleNameArray(jtsTL)
  
  fourNewCols <- str_split(taxaArray,":")
  colData <- data.frame(Sample=taxaArray, TasselIndex = 0:(length(taxaArray)-1L), row.names = taxaArray,
                        matrix(unlist(fourNewCols), nrow = length(fourNewCols), byrow=T))
}


## Constructor for GRanges (GenomicRanges) class object
genomicRanges <- function(genoTable) {
   if(is(genoTable,"GenotypeTable")) {
        jtsPL <- positions(genoTable)@jtsPositionList
    } else if(genoTable %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
      jtsPL <- genoTable$positions()
    } else if(genoTable %instanceof% "net.maizegenetics.dna.map.PositionList") {
      jtsPL <- genoTable
    } else {
        stop("Object is not of \"GenotypeTable\" class")
    }
    
    genoPositionVector <- J("net/maizegenetics/plugindef/GenerateRCode")$genotypeTableToPositionListOfArrays(jtsPL)
    
  
    gr2 <- GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(genoPositionVector$chromosomes),
      ranges = IRanges::IRanges(start = genoPositionVector$startPos, end = genoPositionVector$startPos),
      strand = S4Vectors::Rle(genoPositionVector$strand),
      tasselIndex = 0:(length(genoPositionVector$altAllele)-1L),
      refAllele = genoPositionVector$refAllele,
      altAllele = genoPositionVector$altAllele
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


## Create Summarized Experiment from a TASSEL Genotype Table
snpMatrixFromGenotypeTable <- function(genotypeTable) {
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
  aMatrix <- matrix(genoCallIntArray,length(genomicRangesDF))
  colnames(aMatrix) <- sampleDF[["Sample"]]
  rownames(aMatrix) <- paste0(seqnames(genomicRangesDF),":",ranges(genomicRangesDF))
  as(aMatrix, "SnpMatrix")
}




