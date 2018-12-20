#---------------------------------------------------------------------
# Script Name:   AllClasses.R
# Description:   All classes for rTASSEL
# Author:        Brandon Monier & Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-03 at 17:58:46
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for 
#    TASSEL classes
#--------------------------------------------------------------------





#other methods
#taxa -> vector
#phenotype -> dataframe or tassel obj in wrapper
#genotype -> dataframe or tassel obj in wrapper
#position -> granges or or tassel obj in wrapper



#--------------------------------------------------------------------
# TasselGenotypePhenotype Class
#--------------------------------------------------------------------
setClass(
    Class = "TasselGenotypePhenotype",
    representation = representation(
        name = "character",
        jTasselObj = "jobjRef",
        jTaxaList = "jobjRef",
        jPositionList = "jobjRef",
        jGenotypeTable = "jobjRef",
        jPhenotypeTable = "jobjRef"
    )
)

setMethod(
    f = "show",
    signature = "TasselGenotypePhenotype",
    definition = function(object) {
        cat("A TasselGenotypePhenotype Dataset\n")
        cat("  Class..............", object@name, "\n")
        if (!is.jnull(object@jTaxaList)) {
            cat("  Taxa...............", as.character(object@jTaxaList$size()), "\n")
        } else {
            cat("  Taxa...............", "NA", "\n")
        }
        if (!is.jnull(object@jPositionList)) {
            cat("  Positions..........", as.character(object@jPositionList$numberOfSites()), "\n")
        } else {
            cat("  Positions..........", "NA", "\n")
        }
        if (!is.jnull(object@jTaxaList) & !is.jnull(object@jPositionList)) {
            cat("  Taxa x Positions...", as.numeric(object@jTaxaList$size()) * as.numeric(object@jPositionList$numberOfSites()), "\n")
        } else {
            cat("  Taxa x Positions...", "NA", "\n")
        }
        cat("---\n")
        if (!is.jnull(object@jGenotypeTable)) {
            cat("  Genotype Table..... [x]\n")
        } else {
            cat("  Genotype Table..... [ ]\n")
        }
        if (!is.jnull(object@jPhenotypeTable)) {
            cat("  Phenotype Table.... [x] Traits:", object@jPhenotypeTable$getTableColumnNames(),"\n")
        } else {
            cat("  Phenotype Table.... [ ]\n")
        }
    }
)





# "get" functions

## Get Taxa
getTaxaList <- function(jtsObject) {
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

## Get Positions
getPositionList <- function(jtsObject) {
  if(jtsObject %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
    return(jtsObject$positions())
  } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
    return(jtsObject$genotypeTable()$positions())
  } else {
    return(rJava::.jnull())
  }
}

## Get a GenotypeTable
getGenotypeTable <- function(jtsObject) {
  if(is(jtsObject, "TasselGenotypePhenotype")) {
    return(jtsObject@jGenotypeTable)
  }
  if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
  if(jtsObject %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
    return(jtsObject)
  } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
    return(jtsObject$genotypeTable())
  } else {
    return(rJava::.jnull())
  }
}

## Get a Phenotype object
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



# TasselGenotypePhenotype Object constructors

## main constructor 
.tasselObjectConstructor <- function(jTasselObj) {
  new(
    Class = "TasselGenotypePhenotype",
    name = "TasselGenotypePhenotype",
    jTasselObj = jTasselObj,
    jTaxaList = getTaxaList(jTasselObj),
    jPositionList = getPositionList(jTasselObj),
    jGenotypeTable = getGenotypeTable(jTasselObj),
    jPhenotypeTable = getPhenotypeTable(jTasselObj)
  )
}

## Constructor for GenotypeTable class object
readGenotypeTable <- function(path) {
  .tasselObjectConstructor(rJava::.jcall(
      "net/maizegenetics/dna/snp/ImportUtils",
      "Lnet/maizegenetics/dna/snp/GenotypeTable;",
      "read",
      path
    )
  )
}

## Constructor for PhenotypeTable class object
readPhenotypeTable <- function(path) {
  jObj <- new(J("net.maizegenetics.phenotype.PhenotypeBuilder"))$fromFile(path)
  .tasselObjectConstructor(jObj$build()$get(0L))
}


## TODO Ed
#' @title Constructor for GenotypePhenotype combined object

#' @description Creates a Java GenotypePhenotype object which is used for 
#'    \code{TasselGenotypePhenotype} object construction
#' 
#' @param genoPath a path to a genotype file (e.g. VCF, hmp, etc.) or TASSEL Genotype Obj
#' @param phenoDF a path, a data frame of phenotypic data, or TASSEL Phenotype Obj

# Perhaps add an if statement to determine if phenotype is a data frame or a
#   path
readGenotypePhenotype <- function(genoPathOrObj, phenoPathDFOrObj) {
    genoObj <- getGenotypeTable(genoPathOrObj)
    if(is.jnull(genoObj)) {
      genoObj <- getGenotypeTable(readGenotypeTable(genoPath))
    }
    phenoObj <- getPhenotypeTable(phenoPathDFOrObj)
    if(is.jnull(phenoObj) & is.data.frame(phenoPathDFOrObj)) {
      phenoObj <- createTasselPhenotypeFromDataFrame(phenoPathDFOrObj)
    } else {
      phenoObj <- new(J("net/maizegenetics/phenotype/PhenotypeBuilder"))$fromFile(phenotypePath)$build()$get(0L)
    }
      
    
    t <- .tasselObjectConstructor(
        new(J("net.maizegenetics.phenotype.GenotypePhenotypeBuilder"))
        $genotype(genoObj)$phenotype(phenoObj)$intersect()$build()
    )
}





