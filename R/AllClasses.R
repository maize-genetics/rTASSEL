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
    cat("Taxa: ",object@jtsGenotypeTable$numberOfSites(), " Sites: ",object@jtsGenotypeTable$numberOfTaxa(),"\n")
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


