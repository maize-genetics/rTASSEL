#---------------------------------------------------------------------
# Script Name:   AllClasses.R
# Description:   All classes for rTASSEL
# Author:        Brandon Monier & Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-03 at 17:58:46
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for 
#    TASSEL classes
#--------------------------------------------------------------------

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

## Get positions for GenotypeTable class objects
setMethod(
  f = "positions",
  signature = "GenotypeTable",
  definition = function(object) {
    new("PositionList",name="TASSEL Position List", jtsPositionList=object@jtsGenotypeTable$positions())
  }
)

## Get taxa for GenotypeTable class objects
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

## Show positions lists
setMethod(
  f = "show",
  signature = "PositionList",
  definition = function(object) {
    cat("Position List Name: ",object@name,"\n")
    cat(is(object)," wraps ", show(object@jtsPositionList) ,"\n")
    cat("Sites: ",object@jtsPositionList$size(),"\n")
  }
)

## A R Wrapper for the TaxaList class
setClass(
  Class = "TaxaList",
  representation = representation(
    name = "character",
    jtsTaxaList = "jobjRef"
  )
  #todo - this class seems like it should inherit from jobjRef
  #contains = "jobjRef"
)

## Show method for TaxaList class objects
setMethod(
  f = "show",
  signature = "TaxaList",
  definition = function(object) {
    cat("Taxa List Name: ",object@name,"\n")
    cat(is(object)," wraps ", show(object@jtsTaxaList) ,"\n")
    cat("Taxa: ",object@jtsTaxaList$size(),"\n")
  }
)
