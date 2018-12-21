#---------------------------------------------------------------------
# Script Name:   PhenotypeFunctions.R
# Description:   Functions to create TASSEL Phenotype
# Author:        Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-21 at 15:17:07
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for
#    TASSEL classes
#--------------------------------------------------------------------


## Get a Phenotype object - not exported
getPhenotypeTable <- function(jtsObject) {
  if(is(jtsObject, "TasselGenotypePhenotype")) {
    return(jtsObject@jPhenotypeTable)
  }
  if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
  if(jtsObject %instanceof% "net.maizegenetics.phenotype.Phenotype") {
    return(jtsObject)
  } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
    return(jtsObject$phenotype())
  } else {
    return(rJava::.jnull())
  }
}

#' @title Wrapper function of TasselGenotypePhenotype class for phenotype
#'    data
#' 
#' @description This function is a wrapper for the 
#'    \code{TasselGenotypePhenotype} class. It is used for storing phenotype
#'    information into a class object.
#' 
#' @name readPhenotypeTable
#' @rdname readPhenotypeTable
#' 
#' @param path a phenotype data path or \code{R} data frame
#' 
#' @importFrom rJava J
#' @export
readPhenotypeTable <- function(path) {
  jObj <- new(J("net.maizegenetics.phenotype.PhenotypeBuilder"))$fromFile(path)
  .tasselObjectConstructor(jObj$build()$get(0L))
}


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


