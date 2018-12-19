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
        cat("Class:   ", object@name, "\n")
        if (!is.jnull(object@jTaxaList)) {
            cat("Taxa:    ", as.character(object@jTaxaList$size()), "\n")
        } else {
            cat("Taxa:    ", "NA", "\n")
        }
        if (!is.jnull(object@jPositionList)) {
            cat("Positions:", as.character(object@jPositionList$numberOfSites()), "\n")
        } else {
            cat("Postions: NA", "\n")
        }
    }
)
# ## Show method for TasselObjWrapper class objects
# setMethod(
#   f = "show",
#   signature = "TasselGenotypePhenotype",
#   definition = function(object) {
#     cat("GenotypePhenotype Name: ",object@name,is(object)," wraps ", show(object@jtsObj),"\n")
#     
#     if(!is.null(getTaxaList(object))) {
#       cat("Taxa: ",object@jtsObj$size(),"\n")
#     }
#     #repeat about for position
#     cat("Sites: ",object@jtsPositionList$size(),"\n")
#     #repeat about for call table
#     #repeat for phenotypes
#   }
# )

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

#can return null if missing, otherwise
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

getPositionList <- function(jtsObject) {
  if(jtsObject %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
    return(jtsObject$positions())
  } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
    return(jtsObject$genotypeTable()$positions())
  } else {
    return(rJava::.jnull())
  }
}

getGenotypeTable <- function(jtsObject) {
  if(jtsObject %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
    return(jtsObject)
  } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
    return(jtsObject$genotypeTable())
  } else {
    return(rJava::.jnull())
  }
}

getPhenotypeTable <- function(jtsObject) {
  if(jtsObject %instanceof% "net.maizegenetics.phenotype.Phenotype") {
    return(jtsObject)
  } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
    return(jtsObject$genotypeTable())
  } else {
    return(rJava::.jnull())
  }
}

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

## Constructor for GenotypeTable class object - Terry
readPhenotypeTable <- function(path) {
  jObj <- new(J("net.maizegenetics.phenotype.PhenotypeBuilder"))
  jObj <- jObj$fromFile(path)
  .tasselObjectConstructor(jObj$build()$get(0L))
}

## Constructor for GenotypePhenotype
# TODO - Brandon
readGenotypePhenotype <- function(genoPath, phenoPath) {
    .tasselObjectConstructor(
        new(J("net.maizegenetics.phenotype.GenotypePhenotypeBuilder"))$genotype(genoPath)$phenotype(phenoPath)$intersect()$build()
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

