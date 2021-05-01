#---------------------------------------------------------------------
# Script Name:   AllClasses.R
# Description:   All classes for rTASSEL
# Author:        Brandon Monier & Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2019-09-04 at 14:13:28
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript produce wrapper classes for
#    TASSEL class objects in Java
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# TasselGenotypePhenotype class and constructors
#--------------------------------------------------------------------

#' @title TasselGenotypePhenotype Class
#'
#' @description Class \code{TasselGenotypePhenotype} defines a \code{rTASSEL}
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


#' @title Show method TasselGenotypePhenotype objects
#'
#' @description Prints out information related taxa, positions, genotype, and
#'    phenotype information.
#'
#' @param object a \code{TasselGenotypePhenotype} class object
#'
#' @rdname TasselGenotypePhenotype-class
#' @aliases show,TasselGenotypePhenotype-method
#'
#' @importFrom rJava .jnull
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
            taxaIDs <- object@jPhenotypeTable$getTableColumnNames()
            taxaCutOff <- 5
            if (length(taxaIDs) > taxaCutOff) {
                taxaRem <- length(taxaIDs) - taxaCutOff
                taxaIDs <- taxaIDs[seq_len(taxaCutOff)]
                remMsg <- paste0("... with ", taxaRem, " more IDs\n")
                cat("  Phenotype Table.... [x]\n")
                cat("---\n")
                cat("  Traits:", taxaIDs, remMsg)
            } else {
                cat("  Phenotype Table.... [x]\n")
                cat("---\n")
                cat("  Traits:", taxaIDs, "\n")
            }
        } else {
            cat("  Phenotype Table.... [ ]\n")
        }
    }
)


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


## get TASSEL class - not exported (house keeping)
#' @importFrom methods is
#' @importFrom rJava .jstrVal
.getTASSELClass <- function(object, tasselClassName, throwErrorOnNull = TRUE) {
    jtsObject <- switch(
        tasselClassName,
        "GenotypePhenotype" = getGenotypePhenotype(object),
        "GenotypeTable" = getGenotypeTable(object),
        "Phenotype" = getPhenotypeTable(object),
        "TaxaList" = getTaxaList(object),
        "PositionList" = getPositionList(object)
    )
    if(is.null(jtsObject)) {
        stop("Unknown TASSEL class:", tasselClassName)
    }
    if(throwErrorOnNull & is.jnull(jtsObject)) {
        errObj <- if(is(object,'jobjRef')) {
            rJava::.jstrVal(object)
        } else {
            class(object)
        }
        stop(errObj," does not contain a TASSEL ",tasselClassName," object")
    }
    jtsObject
}


#--------------------------------------------------------------------
# FactorTable class and constructors
#--------------------------------------------------------------------

#' @title FactorTable Class
#'
#' @description Class \code{FactorTable} defines a \code{rTASSEL}
#'    Class for storing TASSEL 6 \code{FactorTable} objects.
#'
#' @name FactorTable-class
#' @rdname FactorTable-class
#' @exportClass TasselGenotypePhenotype
setClass(
    Class = "FactorTable",
    representation = representation(
        jFactorTable = "jobjRef"
    )
)


#' @title Show method FactorTable objects
#'
#' @description Prints out information related taxa, positions, genotype, and
#'    phenotype information.
#'
#' @param object a \code{FactorTable} class object
#'
#' @rdname FactorTable-class
#' @aliases show,FactorTable-method
#'
#' @importFrom rJava .jnull
setMethod(
    f = "show",
    signature = "FactorTable",
    definition = function(object) {
        cat("A FactorTable Dataset\n")
        cat("  Class..........:", class(object), "\n")
        cat("  Taxa...........:", object@jFactorTable$numTaxa(), "\n")
        cat("  Ref. Ranges....:", object@jFactorTable$numFactors(), "\n")
    }
)


