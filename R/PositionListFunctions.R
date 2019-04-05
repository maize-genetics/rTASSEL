#---------------------------------------------------------------------
# Script Name:   PositionListFunctions.R
# Description:   Support working with TASSEL PositionList and GenomicRanges
# Author:        Brandon Monier & Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2019-04-04 at 22:50:01
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions
#    necessary for extracting position list data for
#    `SummarizedExperiment` objects.
#--------------------------------------------------------------------

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
