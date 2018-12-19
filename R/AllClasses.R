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

## Make a GenotypeTable class
setClass(
  Class = "GenotypeTable",
  representation = representation(
    name = "character",
    jtsGenotypeTable = "jobjRef"
  )
  #todo - this class seems like it should inherit from jobjRef
  #contains = "jobjRef"
)

## Display overview when object is called
setMethod(
  f = "show",
  signature = "GenotypeTable",
  definition = function(object) {
    cat("Genotype Table Name: ",object@name,"\n")
    cat(is(object)," wraps ", show(object@jtsGenotypeTable) ,"\n")
    cat("Sites: ",object@jtsGenotypeTable$numberOfSites(), " Taxa: ",object@jtsGenotypeTable$numberOfTaxa(),"\n")
  }
)

## Get positions for GenotypeTable class objects
setMethod(
  f = "positions",
  signature = "GenotypeTable",
  definition = function(object) {
    new("PositionList",name="TASSEL Position List", jtsPositionList=object@jtsGenotypeTable$positions())
  }
)

## Get taxa for GenotypeTable class objects
setMethod(
  f = "taxa",
  signature = "GenotypeTable",
  definition = function(object) {
    new("TaxaList",name="TASSEL Taxa List", jtsTaxaList=object@jtsGenotypeTable$taxa())
  }
)

## Constructor for GenotypeTable class object
readGenotypeTable <- function(path) {
  new(
    Class = "GenotypeTable",
    #todo split path to only grab the filename
    name = paste0("PATH:",path),
    jtsGenotypeTable = rJava::.jcall(
      "net/maizegenetics/dna/snp/ImportUtils",
      "Lnet/maizegenetics/dna/snp/GenotypeTable;",
      "read",
      path
    )
  )
}

## A R Wrapper for the PositionList class
setClass(
  Class = "PositionList",
  representation = representation(
    name = "character",
    jtsPositionList = "jobjRef"
  )
  #todo - this class seems like it should inherit from jobjRef
  #contains = "jobjRef"
)

## Show positions lists
setMethod(
  f = "show",
  signature = "PositionList",
  definition = function(object) {
    cat("Position List Name: ",object@name,"\n")
    cat(is(object)," wraps ", show(object@jtsPositionList) ,"\n")
    cat("Sites: ",object@jtsPositionList$size(),"\n")
  }
)

## A R Wrapper for the TaxaList class
setClass(
  Class = "TaxaList",
  representation = representation(
    name = "character",
    jtsTaxaList = "jobjRef"
  )
  #todo - this class seems like it should inherit from jobjRef
  #contains = "jobjRef"
)

## Show method for TaxaList class objects
setMethod(
  f = "show",
  signature = "TaxaList",
  definition = function(object) {
    cat("Taxa List Name: ",object@name,"\n")
    cat(is(object)," wraps ", show(object@jtsTaxaList) ,"\n")
    cat("Taxa: ",object@jtsTaxaList$size(),"\n")
  }
)

## A R Wrapper for the TasselObjWrapper class
setClass(
  Class = "TasselObjWrapper",
  representation = representation(
    name = "character",
    jtsObj = "jobjRef"
  )
  #todo - this class seems like it should inherit from jobjRef
  #contains = "jobjRef"
)

## Show method for TasselObjWrapper class objects
setMethod(
  f = "show",
  signature = "TasselObjWrapper",
  definition = function(object) {
    #SummarizeTaxa
    
    #Summarize Genotypic side if applicable
    cat("GenotypePhenotype Name: ",object@name,"\n")
    cat(is(object)," wraps ", show(object@jtsObj) ,"\n")
    cat("Taxa: ",object@jtsObj$size(),"\n")
    #Summarize Phenotypic side if applicable
  }
)

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
        jtsGenotypePhenotypeTable = "jobjRef"
    )
)

#--------------------------------------------------------------------
# PSEUDO CODE - TasselGenotypePhenotype Object and constructors
#--------------------------------------------------------------------
TasselGenotypePhenotype <- function(genotype, phenotype) {
    if (missing(genotype) & missing(phenotype)) {
        stop("Need at least one path.")
    } else if (missing(genotype)) {
        if (isPathOrRObj(phenotype) == "path") {
            javaObj <- "createJPhenotypeFromPath(phenotype)"
        } else {
            javaObj <- "createJPhenotypeFromRObj(phenotype)"
        }
    } else if (missing(phenotype)) {
        if (isPathOrRObj(genotype) == "path") {
            javaObj <- "createJGenotypeFromPath(genotype)"
        } else {
            javaObj <- "createJGenotypeFromRObj(genotype)"
        }
    } else {
        if (isCombinedPathOrRObj(genotype, phenotype) == "path") {
            javaObj <- "createJGenotypePhenotypeFromPath(genotype, phenotype)"
        } else {
            javaObj <- "createJGenotypePhenotypeFromRObj(genotype, phenotype)"
        }
    }
    new(
        Class = "TasselGenotypePhenotype",
        name = "TasselGenotypePhenotype",
        jtsGenotypePhenotypeTable = javaObj
    )
}

# ## TasselGenotypePhenotype Show Method
# setMethod(
#     f = "show",
#     signature = "TasselGenotypePhenotype",
#     definition = function(object) {
#         cat("Class: ", class(object),"\n")
#         cat("\n\nGenotype Info:\n")
#         print(object@genotypeTable)
#         cat("\n\nPhenotype Info:\n")
#         print(object@phenotypeTable)
#     }
# )

# ## TasselGenotypePhenotype Getters - genotypeTable
# setMethod(
#     f = "genotypeTable",
#     signature = "TasselGenotypePhenotype",
#     definition = function(object) {
#         object@genotypeTable
#     }
# )

# ## TasselGenotypePhenotype Getters - phenotypeTable
# setMethod(
#     f = "phenotypeTable",
#     signature = "TasselGenotypePhenotype",
#     definition = function(object) {
#         object@phenotypeTable
#     }
# )

