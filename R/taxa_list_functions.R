## ----
# Get Taxa - not exported (house keeping)
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
#' @rdname taxaList
#' @aliases taxaList,TasselGenotypePhenotype-method
#' @export
setMethod("taxaList", "TasselGenotypePhenotype", function(tasObj) {
    return(getTaxaIDs(tasObj))
})



## Methods for pulling Taxa or Samples - not exported (house keeping)
#' @importFrom rJava J
getTaxaIDs <- function(tasObj) {
    if (!inherits(tasObj, "TasselGenotypePhenotype")) {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jtsTL <- getTaxaList(tasObj)
    rJava::J("net/maizegenetics/plugindef/GenerateRCode")$
        genotypeTableToSampleNameArray(jtsTL)
}


## Get sample ID data frame - not exported (house keeping)
sampleDataFrame <- function(tasObj) {
    taxaArray <- getTaxaIDs(tasObj)

    tibble::tibble(
        Sample = taxaArray,
        TasselIndex = 0:(length(taxaArray) - 1L)
    )
}


