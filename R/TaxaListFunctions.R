#---------------------------------------------------------------------
# Script Name:   TaxaListFunctions.R
# Description:   Functions to support TaxaList or Samples
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2019-04-04 at 16:32:24
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#   The main purpose of this Rscript is to house functions
#   necessary for extracting TASSEL taxa lists from TASSEL objects.
#--------------------------------------------------------------------

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


## Methods for pulling Taxa or Samples - not exported (house keeping)
sampleVectorFromTassel <- function(tasObj) {
    jtsTL <- getTaxaList(tasObj)
    rJava::J("net/maizegenetics/plugindef/GenerateRCode")$
        genotypeTableToSampleNameArray(jtsTL)
}


## Get sample ID data frame - not exported (house keeping)
sampleDataFrame <- function(tasObj) {
    taxaArray <- sampleVectorFromTassel(tasObj)
    # fourNewCols <- stringr::str_split(taxaArray, ":")
    colData <- data.frame(
        row.names = taxaArray,
        Sample = taxaArray,
        TasselIndex = 0:(length(taxaArray) - 1L)
        # DEBUG = matrix(
        #     unlist(fourNewCols),
        #     nrow = length(fourNewCols),
        #     byrow = TRUE
        # )
    )
    return(colData)
}
