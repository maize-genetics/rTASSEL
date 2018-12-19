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

## Show method for TasselObjWrapper class objects
setMethod(
  f = "show",
  signature = "TasselGenotypePhenotype",
  definition = function(object) {
    cat("GenotypePhenotype Name: ",object@name,is(object)," wraps ", show(object@jtsObj),"\n")
    
    if(!is.null(getTaxaList(object))) {
      cat("Taxa: ",object@jtsObj$size(),"\n")
    }
    #repeat about for position
    cat("Sites: ",object@jtsPositionList$size(),"\n")
    #repeat about for call table
    #repeat for phenotypes
  }
)

.tasselObjectConstructor <- function(javaObj) {
  new(
    Class = "TasselGenotypePhenotype",
    name = "TasselGenotypePhenotype",
    jTasselObj = javaObj,
    jTaxaList = getTaxaList(jTasselObj),
    # continue for rest
  )
}

#can return null if missing, otherwise 
getTaxaList <- function(jtsObject) {
  if(jtsObject %instanceof% TaxaList) {
    return jtsObject
  }
  if(jtsObject %instanceof% GenotypeTable) {
    return jtsObject$taxa()
  }
  if(jtsObject %instanceof% Phenotype) {
    return jtsObject$taxa()
  }
  #something for GenotypePhenotype
  if(jtsObject %instanceof% GenotypePhenotype) {
    return jtsObject$genotypeTable()$taxa()
  }
  return NULL
}

getPositionList <- function(jtsObject) {
  if(jtsObject %instanceof% GenotypeTable) {
    return jtsObject$positions()
  }
  if(jtsObject %instanceof% GenotypePhenotype) {
    return jtsObject$genotypeTable()$positions()
  }
  return NULL
}

getGenotypeTable <- function(jtsObject) {
  if(jtsObject %instanceof% GenotypeTable) {
    return jtsObject
  }
  if(jtsObject %instanceof% GenotypePhenotype) {
    return jtsObject$genotypeTable()
  }
  return NULL
}

getPhenotypeTable <- function(jtsObject) {
  if(jtsObject %instanceof% Phenotype) {
    return jtsObject
  }
  #something for GenotypePhenotype
  if(jtsObject %instanceof% GenotypePhenotype) {
    return jtsObject$genotypeTable()
  }
  return NULL
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
  .tasselObject(rJava::.jcall(
      "net/maizegenetics/dna/snp/ImportUtils",
      "Lnet/maizegenetics/dna/snp/GenotypeTable;",
      "read",
      path
    )
  )
}

## Constructor for GenotypeTable class object
readPhenotypeTable <- function(path) {
  # .tasselObject(rJava::.jcall(
  #     "net/maizegenetics/dna/snp/ImportUtils",
  #     "Lnet/maizegenetics/dna/snp/GenotypeTable;",
  #     "read",
  #     path
  #   )
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

