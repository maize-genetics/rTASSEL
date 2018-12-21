#---------------------------------------------------------------------
# Script Name:   AllClasses.R
# Description:   All classes for rTASSEL
# Author:        Brandon Monier & Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2018-12-21 at 16:11:58
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for 
#    TASSEL classes
#--------------------------------------------------------------------

# TODO All
#   other methods
#   taxa -> vector
#   phenotype -> dataframe or tassel obj in wrapper
#   genotype -> dataframe or tassel obj in wrapper
#   position -> granges or or tassel obj in wrapper



#--------------------------------------------------------------------
# TasselGenotypePhenotype class and constructors
#--------------------------------------------------------------------

#' @title TasselGenotypePhenotype Class
#' 
#' @description Class \code{TasselGenotypePhenoty} defines a \code{rTASSEL}
#'    Class for storing TASSEL genotype and phenotype objects.
#'
#' @name TasselGenotypePhenotype-class
#' @rdname TasselGenotypePhenotype-class
#' @exportClass TasselGenotypePhenotype
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


#' @title Wrapper function of TasselGenotypePhenotype class for genotype
#'    data
#' 
#' @description This function is a wrapper for the 
#'    \code{TasselGenotypePhenotype} class. It is used for storing genotype
#'    information into a class object.
#' 
#' @name readGenotypeTable
#' @rdname readGenotypeTable
#' 
#' @param path a genotype data path (e.g. \code{*.VCF, *.hmp}, etc.)
#' 
#' @export
readGenotypeTable <- function(path) {
  .tasselObjectConstructor(rJava::.jcall(
      "net/maizegenetics/dna/snp/ImportUtils",
      "Lnet/maizegenetics/dna/snp/GenotypeTable;",
      "read",
      path
    )
  )
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
#' @export
readPhenotypeTable <- function(path) {
  jObj <- new(J("net.maizegenetics.phenotype.PhenotypeBuilder"))$fromFile(path)
  .tasselObjectConstructor(jObj$build()$get(0L))
}


#' @title Wrapper function of TasselGenotypePhenotype class for 
#'    GenotypePhenotype combined data

#' @description Creates a Java GenotypePhenotype object which is used for 
#'    \code{TasselGenotypePhenotype} object construction
#' 
#' @name readGenotypePhenotype
#' @rdname readGenotypePhenotype
#' 
#' @param genoPathOrObj a path to a genotype file (e.g. VCF, hmp, etc.) or 
#'    TASSEL Genotype Obj
#' @param phenoPathDFOrObj a path, a data frame of phenotypic data, or TASSEL 
#'    Phenotype Obj
#'
#' @export
readGenotypePhenotype <- function(genoPathOrObj, phenoPathDFOrObj) {
    genoObj <- getGenotypeTable(genoPathOrObj)
    if(is.jnull(genoObj)) {
      genoObj <- getGenotypeTable(readGenotypeTable(genoPathOrObj))
    }
    phenoObj <- getPhenotypeTable(phenoPathDFOrObj)
    if(is.jnull(phenoObj) & is.data.frame(phenoPathDFOrObj)) {
      phenoObj <- createTasselPhenotypeFromDataFrame(phenoPathDFOrObj)
    } else {
      phenoObj <- new(J("net/maizegenetics/phenotype/PhenotypeBuilder"))$fromFile(phenoPathDFOrObj)$build()$get(0L)
    }

    t <- .tasselObjectConstructor(
        new(J("net.maizegenetics.phenotype.GenotypePhenotypeBuilder"))
        $genotype(genoObj)$phenotype(phenoObj)$intersect()$build()
    )
}


#' @title Show method TasselGenotypePhenotype objects
#' 
#' @description Prints out information related taxa, positions, genotype, and
#'    phenotype information.
#'
#' @rdname show-methods
#' @aliases show,TasselGenotypePhenotype-method
#' @exportMethod show
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



#--------------------------------------------------------------------
# "get" functions for Core functions
#--------------------------------------------------------------------

## Get Taxa - not exported
getGenotypePhenotype <- function(jtsObject) {
  if(is(jtsObject, "TasselGenotypePhenotype")) {
    jtsObject <- jtsObject@jTasselObj
  }
  if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
  if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
    return(jtsObject)
  } else {
    return(rJava::.jnull())
  }
}

## Get Taxa - not exported
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

## Get a GenotypeTable - not exported
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



#--------------------------------------------------------------------
# Core functions for TasselGenotypePhenotype class objects
#--------------------------------------------------------------------

## main constructor for TasselGenotypePhenotype objects - not exported
.tasselObjectConstructor <- function(jTasselObj) {
  tobj <- new(
    Class = "TasselGenotypePhenotype",
    name = "TasselGenotypePhenotype",
    jTasselObj = jTasselObj,
    jTaxaList = getTaxaList(jTasselObj),
    jPositionList = getPositionList(jTasselObj),
    jGenotypeTable = getGenotypeTable(jTasselObj),
    jPhenotypeTable = getPhenotypeTable(jTasselObj)
  )
  if(is.jnull(tobj@jTaxaList) & is.jnull(tobj@jPositionList) & 
     is.jnull(tobj@jGenotypeTable) & is.jnull(tobj@jPhenotypeTable)) {
    return (NULL)
  }
  tobj
}


## get TASSEL class - not exported
.getTASSELClass <- function(object, tasselClassName, throwErrorOnNull = TRUE) {
  jtsObject <- switch(tasselClassName,
                      "GenotypePhenotype" = getGenotypePhenotype(object),
                      "GenotypeTable" = getGenotypeTable(object),
                      "Phenotype" = getPhenotypeTable(object),
                      "TaxaList" = getTaxaList(object),
                      "PositionList" = getPositionList(object)
                      )
  if(is.null(jtsObject)) {
    stop("Unknown TASSEL class:",tasselClassName)
  }
  if(throwErrorOnNull & is.jnull(jtsObject)) {
    errObj <- if(is(object,'jobjRef')) .jstrVal(object) else class(object)
    stop(errObj," does not contain a TASSEL ",tasselClassName," object")
  }
  jtsObject
}
