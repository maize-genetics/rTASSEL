#---------------------------------------------------------------------
# Script Name:   PositionListFunctions.R
# Description:   Support working with TASSEL PositionList and GenomicRanges
# Author:        Brandon Monier & Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-21 at 15:16:56
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for
#    TASSEL classes
#--------------------------------------------------------------------



## Get Positions - not exported
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

## Constructor for GRanges (GenomicRanges) class object
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
