##--------------------------------------------------------------------
# Script Name:   rjava_tests.R
# Description:   Various tests with rJava
# Author:        Brandon Monier & Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-03 at 17:58:46
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for 
#    TASSEL classes
#--------------------------------------------------------------------

# Preamble 

## Make a GenotypeTable class
setClass(
  Class = "GenotypeTable",
  representation = representation(
    name = "character",
    jtsGenotypeTable = "jobjRef"
  )
  #todo - this class seems like it should inherit from jobjRef
  #contains = "jobjRef"
)

## Display overview when object is called
setMethod(
  f = "show",
  signature = "GenotypeTable",
  definition = function(object) {
    cat("Genotype Table Name: ",object@name,"\n")
    cat(is(object)," wraps ", show(object@jtsGenotypeTable) ,"\n")
    cat("Sites: ",object@jtsGenotypeTable$numberOfSites(), " Taxa: ",object@jtsGenotypeTable$numberOfTaxa(),"\n")
  }
)

setMethod(
  f = "positions",
  signature = "GenotypeTable",
  definition = function(object) {
    new("PositionList",name="TASSEL Position List", jtsPositionList=object@jtsGenotypeTable$positions())
  }
)

setMethod(
  f = "taxa",
  signature = "GenotypeTable",
  definition = function(object) {
    new("TaxaList",name="TASSEL Taxa List", jtsTaxaList=object@jtsGenotypeTable$taxa())
  }
)

## Constructor for GenotypeTable class object
readGenotypeTable <- function(path) {
  new(
    Class = "GenotypeTable",
    #todo split path to only grab the filename
    name = paste0("PATH:",path),
    jtsGenotypeTable = rJava::.jcall(
      "net/maizegenetics/dna/snp/ImportUtils",
      "Lnet/maizegenetics/dna/snp/GenotypeTable;",
      "read",
      path
    )
  )
}

## A R Wrapper for the PositionList class
setClass(
  Class = "PositionList",
  representation = representation(
    name = "character",
    jtsPositionList = "jobjRef"
  )
  #todo - this class seems like it should inherit from jobjRef
  #contains = "jobjRef"
)

setMethod(
  f = "show",
  signature = "PositionList",
  definition = function(object) {
    cat("Position List Name: ",object@name,"\n")
    cat(is(object)," wraps ", show(object@jtsPositionList) ,"\n")
    cat("Sites: ",object@jtsPositionList$size(),"\n")
  }
)

setClass(
  Class = "TaxaList",
  representation = representation(
    name = "character",
    jtsTaxaList = "jobjRef"
  )
  #todo - this class seems like it should inherit from jobjRef
  #contains = "jobjRef"
)

setMethod(
  f = "show",
  signature = "TaxaList",
  definition = function(object) {
    cat("Taxa List Name: ",object@name,"\n")
    cat(is(object)," wraps ", show(object@jtsTaxaList) ,"\n")
    cat("Taxa: ",object@jtsTaxaList$size(),"\n")
  }
)

## Constructor for GenotypeTable class object
sampleDataFrame <- function(jtsGenoTableOrTaxaList) {
  if(is(jtsGenoTableOrTaxaList,"GenotypeTable")) {
    jtsTL <- taxa(jtsGenoTableOrTaxaList)@jtsTaxaList
  } else if(is(jtsGenoTableOrTaxaList,"TaxaList")) {
    jtsTL <- jtsGenoTableOrTaxaList@jtsTaxaList
  } else {
    jtsTL <- jtsGenoTableOrTaxaList
  }
  taxaArray <- c()
  for(i in 1:jtsTL$size()) {
    #why do I have to do -1L
    taxaArray[i] = toString(jtsTL$taxaName(i-1L))
  }
  colData <- data.frame(Sample=taxaArray)
  colData
}


## Constructor for GRanges (GenomicRanges) class object
sampleGenomicRanges <- function(jtsGenoTable) {
    if(is(jtsGenoTable,"GenotypeTable")) {
        jtsGT <- positions(jtsGenoTable)@jtsPositionList
    } else {
        stop("Object is not of \"GenotypeTable\" class")
    }
    
    numSite <- as.numeric(jtsGT$numberOfSites())
    physPos <- jtsGT$physicalPositions()
    
    cat("Extracting chromosome names for each postion...\n")
    cat("...is there a quicker way to get this? (~ Brandon)\n")
    chrName <- lapply(seq_len(numSite), function(pos) {
        jtsGT$chromosomeName(as.integer(pos - 1))
    })
    chrName <- unlist(chrName)
    
    gr2 <- GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(chrName),
        ranges = IRanges::IRanges(start = physPos, end = physPos)
    )
    return(gr2)
}