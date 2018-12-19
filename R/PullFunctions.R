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

createTasselDataSet <- function(...) {
  arguments <- list(...)
  jList <- new(J("java/util/ArrayList"))
  for(javaObj in arguments) {
    #check if they are all TASSEL jobj
    if(is(javaObj, "jobjRef") == FALSE) {
      stop(paste0("Object ",javaObj," is not of class"))
    }
    jList$add(new(J("net/maizegenetics/plugindef/Datum"),"FromR",javaObj,NULL))
  }
  new(J("net/maizegenetics/plugindef/DataSet"),jList,NULL)
}

combineTasselGenotypePhenotype <- function(genotypeTable, phenotype) {
  new(J("net.maizegenetics.phenotype.GenotypePhenotypeBuilder"))$genotype(genotypeTable)$phenotype(phenotype)$intersect()$build()
}

convertTableReportToDataFrame <- function(tableReport) {
  tableReportVectors <- J("net/maizegenetics/plugindef/GenerateRCode")$tableReportToVectors(tableReport)
  colNum <- length(tableReportVectors$columnNames)
  aDF <- data.frame(tableReportVectors$dataVector$get(0L))
  for(i in 2:(colNum)) {
    aDF[[i]] <- tableReportVectors$dataVector$get(i-1L)
  }
  colnames(aDF) <- tableReportVectors$columnNames
  aDF
}

dataSetToListOfDataFrame <- function(jtsDataSet) {
  result <- c()
  for(i in 1:(jtsDataSet$getSize())) {
    name <- jtsDataSet$getData(i-1L)$getName()
    result[[name]] <- convertTableReportToDataFrame(jtsDataSet$getData(i-1L)$getData())
  }
  result
}



