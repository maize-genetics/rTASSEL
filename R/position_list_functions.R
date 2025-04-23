## ----
## Get Positions - not exported (house keeping)
getPositionList <- function(jtsObject) {
    if(is(jtsObject, "TasselGenotypePhenotype")) {
        return(jtsObject@jPositionList)
    }
    if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
    if(jtsObject %instanceof% "net.maizegenetics.dna.map.PositionList") {
        return(jtsObject)
    } else if(jtsObject %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
        return(jtsObject$positions())
    } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
        return(jtsObject$genotypeTable()$positions())
    } else {
        return(rJava::.jnull())
    }
}


## ----
## Constructor for GRanges (GenomicRanges) class object - not exported (<TMP>)
genomicRanges <- function(genoTable) {
    jtsPL <- .getTASSELClass(genoTable, "PositionList")

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


## ----
#' @title Get position list metadata from genotype table
#'
#' @description Returns positional data from a \code{TasselGenotypePhenotype}
#'    object
#'
#' @param tasObj A \code{TasselGenotypePhenotype} object
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava J
#' @importFrom rJava new
#'
#' @export
positionList <- function(tasObj) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    if (rJava::is.jnull(tasObj@jGenotypeTable)) {
        stop("`tasObj` must contain genotype data")
    }

    sites <- rJava::new(
        rJava::J("net.maizegenetics.dna.map.PositionListTableReport"),
        getPositionList(tasObj)
    )

    return(tableReportToDF(sites))
}


