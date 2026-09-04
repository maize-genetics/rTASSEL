## ----
# Get Taxa - not exported (house keeping)
getTaxaList <- function(jtsObject) {
    if(is(jtsObject, "TasselGenotypePhenotype")) {
        return(jtsObject@jTaxaList)
    }
    jtsObject <- .unwrapTasselObject(jtsObject)
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
##
## Any object that carries either genotype or phenotype data has a taxa
## list, so this reaches for the list directly rather than going through
## '.resolveTasselInput()' and its component requirements.
#' @importFrom rJava J
getTaxaIDs <- function(tasObj) {
    jtsTL <- getTaxaList(tasObj)

    if (rJava::is.jnull(jtsTL)) {
        rlang::abort(c(
            "No taxa found in input",
            "x" = sprintf(
                "Got an object of class <%s>",
                paste(class(tasObj), collapse = "/")
            )
        ))
    }

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


