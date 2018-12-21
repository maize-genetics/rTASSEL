#---------------------------------------------------------------------
# Script Name:   PushFunctions.R
# Description:   Functions to create TASSEL objects from rJava
# Author:        Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-03 at 17:58:46
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for
#    TASSEL classes
#--------------------------------------------------------------------

source("R/AllClasses.R")

createTasselPhenotypeFromDataFrame <- function(phenotypeDF, attributeTypes = NULL) {
  taxaNames <- as.vector(phenotypeDF$Taxon)
  colnames <- colnames(phenotypeDF)
  notTaxaCols <- colnames[!colnames %in% c("Taxon")]
  if(is.null(attributeTypes)) {
    atttype <- c(rep("data",length(notTaxaCols)))
  } else {
    atttype <- attributeTypes
  }
  jList <- new(J("java/util/ArrayList"))
  for (col_i in notTaxaCols) {
    jList$add(.jarray(phenotypeDF[[col_i]]))
    
  }
  jc <- J("net/maizegenetics/plugindef/GenerateRCode")$createPhenotypeFromRDataFrameElements(taxaNames,notTaxaCols,atttype,jList)
  .tasselObjectConstructor(jc)
}

