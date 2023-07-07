#---------------------------------------------------------------------
# Script Name:   TaxaListFunctions.R
# Description:   Functions to support TaxaList or Samples
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2020-06-16 at 17:27:42
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#   The main purpose of this Rscript is to house functions
#   necessary for extracting TASSEL taxa lists from TASSEL objects.
#--------------------------------------------------------------------

## ----
## Get Taxa - not exported (house keeping)
getTaxaList <- function(jtsObject) {
    if(is(jtsObject, "TasselGenotypePhenotype")) {
        return(jtsObject@jTaxaList)
    }
    if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
    if(jtsObject %instanceof% "net.maizegenetics.taxa.TaxaList") {
        return(jtsObject)
    } else if(jtsObject %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
        return(jtsObject$taxa())
    } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.Phenotype") {
        return(jtsObject$taxa())
    } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
        return(jtsObject$genotypeTable()$taxa())
    } else {
        return(rJava::.jnull())
    }
}


## ----
#' @title Get list of taxa from TASSEL data
#'
#' @description Returns a list of taxa from a \code{TasselGenotypePhenotype}
#'    object
#'
#' @param tasObj A \code{TasselGenotypePhenotype} object
#'
#' @export
taxaList <- function(tasObj) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    return(getTaxaIDs(tasObj))
}



## Methods for pulling Taxa or Samples - not exported (house keeping)
#' @importFrom rJava J
getTaxaIDs <- function(tasObj) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jtsTL <- getTaxaList(tasObj)
    rJava::J("net/maizegenetics/plugindef/GenerateRCode")$
        genotypeTableToSampleNameArray(jtsTL)
}


## Get sample ID data frame - not exported (house keeping)
sampleDataFrame <- function(tasObj) {
    taxaArray <- getTaxaIDs(tasObj)

    S4Vectors::DataFrame(
        Sample = taxaArray,
        TasselIndex = 0:(length(taxaArray) - 1L)
    )
}


