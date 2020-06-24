#--------------------------------------------------------------------
# Script Name:   AssociationFunctions.R
# Description:   General functions for running association analysis
# Author:        Brandon Monier
# Created:       2019-03-29 at 13:51:02
# Last Modified: 2020-06-24 at 09:50:27
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house all necessary
#    general functions and wrappers for TASSEL association analyses.
#--------------------------------------------------------------------

#' @title R interface for TASSEL's association methods
#'
#' @description This function acts as a front-end for TASSEL's extensive
#'    association analysis methods. Using this function, users can run
#'    the following TASSEL association methods:
#'    \itemize{
#'      \item{best linear unbiased estimates (BLUEs)}
#'      \item{generalized linear model (GLM)}
#'      \item{mixed linear model}
#'      \item{Fast association (Shabalin 2012)}
#'    }
#'
#' @name assocModelFitter
#' @rdname assocModelFitter
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#' @param formula An R-based linear model formula. The general layout of this
#'   formula uses the following TASSEL data scheme:
#'   \code{<data> ~ <factor> and/or <covariate>}. If all traits in a Phenotype
#'   object should be ran, a simplified formula (\code{. ~ .}) can be used.
#'   This scheme can also be used for running all \code{<data>} or
#'   \code{<factor>} and/or \code{<covariate>} data as well. Single variables
#'   are separated witha \code{+} operator. See vignette for further
#'   clarification.
#' @param fitMarkers Should marker data be fitted? If \code{TRUE}, GLM
#'   analysis will be executed. If \code{FALSE}, BLUEs will be calculated.
#'   Defaults to \code{FALSE}.
#' @param kinship Should kinship data be accounted for in the model? If so,
#'   please submit a TASSEL kinship matrix object using the rTASSEL function,
#'   \code{kinshipMatrix}. Defaults to \code{NULL}
#' @param fastAssociation Should TASSEL's Fast Association plugin be used?
#'   Consider setting to \code{TRUE} if you have many phenotypes in your
#'   data set.
#' @param maxP Maximum p-value (0 - 1) to be reported. Currently works with
#'   fast association only. Defaults to a p-value of \code{0.001}
#'   will be used as a threshold. \strong{Note:} p-value parameter will
#'   not be used for BLUE analysis.
#' @param maxThreads Maximum threads to be used when running fast association.
#'   If \code{NULL}, all threads on machine will be used.
#' @param outputFile Output file prefix to be specified in case you want
#'   to write data directly to disk. Highly recommended for large datasets.
#'   If \code{NULL}, no data will be saved to disk. If a character
#' @param minClassSize The minimum acceptable genotype class size. Genotypes
#'   in a class with a smaller size will be set to missing. Defaults to 0.
#' @param biallelicOnly Only test sites that are bi-allelic. The alternative is
#'   to test sites with two or more alleles. Defaults to \code{FALSE}
#' @param appendAddDom If true, additive and dominance effect estimates will
#'   be added to the stats report for bi-allelic sites only. The effect will
#'   only be estimated when the data source is genotype (not a probability).
#'   The additive effect will always be non-negative. Defaults to \code{FALSE}.
#'
#' @return Returns an R list containing \code{DataFrame}-based data frames
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava J
#' @importFrom rJava .jnew
#' @importFrom rJava .jnull
#' @importFrom rJava new
#' @importFrom rlang .data
#' @export
assocModelFitter <- function(tasObj,
                             formula,
                             fitMarkers = FALSE,
                             kinship = NULL,
                             fastAssociation = FALSE,
                             maxP = 0.001,
                             maxThreads = NULL,
                             minClassSize = 0,
                             outputFile = NULL,
                             biallelicOnly = FALSE,
                             appendAddDom = FALSE) {

    # Logic - Check for TasselGenotypePhenotype class
    if (!is(tasObj, "TasselGenotypePhenotype")) {
        stop("tasObj is not of class \"TasselGenotypePhenotype\"")
    }

    # Logic - Check to see if TASSEL object has a phenotype table
    if (rJava::is.jnull(tasObj@jPhenotypeTable)) {
        stop("tasObj does not contain a Phenotype object")
    }

    # Extract formula response and prediction components
    formResp <- all.vars(formula[[2]])
    formPred <- all.vars(formula[[3]])

    # Get all TASSEL object trait metadata
    jtsPheno <- getPhenotypeTable(tasObj)
    phenoAttDf <- extractPhenotypeAttDf(jtsPheno)

    # Add "." variable for whitelisting
    wildCard <- data.frame(
        traitName = ".",
        traitType = "wildcard",
        traitAttribute = "AnyAttribute"
    )
    phenoAttDf <- rbind(phenoAttDf, wildCard)

    # Subset TASSEL object trait types
    tasResp <- subset(
        x = phenoAttDf,
        traitType == "data" |
            traitType == "wildcard"
    )
    tasPred <- subset(
        x = phenoAttDf,
        traitType == "factor" |
            traitType == "covariate" |
            traitType == "wildcard"
    )

    # Logic - Check formula entry
    if (any(!(c(formResp, formPred) %in% phenoAttDf$traitName))) {
        stop("Variables in formula do not match traits in TASSEL object.")
    } else if (any(!(formResp %in% tasResp$traitName))) {
        stop("Only <data> trait types can be implemented as response variables.")
    } else if (any(!(formPred %in% tasPred$traitName))) {
        stop("Only <factor> or <covariate> trait types can be implemented as predictor variables.")
    }

    # Logic - Handle "." variables
    if (all(formResp == ".", formPred == ".")) {
        message("Running all traits...")
        finalResp <- as.vector(tasResp[which(tasResp$traitName != "."), ]$traitName)
        finalPred <- as.vector(tasPred[which(tasPred$traitName != "."), ]$traitName)
    } else if (all(formResp == ".", formPred != ".")) {
        message("Running all <data> traits...")
        finalResp <- as.vector(tasResp[which(tasResp$traitName != "."), ]$traitName)
        finalPred <- formPred
    } else if (all(formResp != ".", formPred == ".")) {
        message("Running all non <data> traits and/or <taxa>...")
        finalResp <- formResp
        finalPred <- as.vector(tasPred[which(tasPred$traitName != "."), ]$traitName)
    } else {
        finalResp <- formResp
        finalPred <- formPred
    }

    # Logic - Handle association analyses
    jRC <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    jTasFilt <- tasPhenoFilter(
        tasObj = tasObj,
        filtObj = c(finalResp, finalPred)
    )
    tmpDF <- as.data.frame(jTasFilt$phenoDf) # check for missing values

    # Logic - Check for out of range p-values
    if (maxP > 1 || maxP < 0) {
        stop("p-value is out of range (0 - 1)")
    }

    # Logic - Convert p-values to Java data types
    if (!is.numeric(maxP)) {
        stop("p-value must be numeric")
    } else if (!fastAssociation) {
        maxP <- as.double(maxP)
    } else {
        maxP <- rJava::.jnew("java/lang/Double", as.character(maxP))
    }

    # Logic - Convert threads to Java data types
    if (!is.null(maxThreads)) {
        maxThreads <- rJava::.jnew("java/lang/Integer", toString(maxThreads))
    } else {
        maxThreads <- rJava::.jnull()
    }

    # Logic - Check output data
    if (!is.null(outputFile)) {
        saveToFile <- TRUE
    } else {
        saveToFile <- FALSE
        outputFile <- rJava::.jnull()
    }

    # Logic - Handle association types and output
    if (!fitMarkers & is.null(kinship)) {
        if (!fastAssociation) {
            if (!rJava::is.jnull(tasObj@jPhenotypeTable)) {
                message("Association Analysis : BLUEs")
                assocOut <- jRC$association(
                    rJava::.jnull(),
                    rJava::.jnull(),
                    jTasFilt$phenotype,
                    rJava::.jnull(),
                    as.integer(minClassSize),
                    biallelicOnly,
                    appendAddDom,
                    saveToFile,
                    outputFile,
                    maxP
                )
                assocType <- "BLUE"
            } else {
                stop("No TASSEL phenotype table was found in TasselGenotypePhenotype object!")
            }
        } else {
            stop("Don't know how to analyze with given parameter inputs.")
        }
    } else if (fitMarkers & is.null(kinship)) {
        if (!fastAssociation) {
            if (!rJava::is.jnull(tasObj@jGenotypeTable)) {
                if (any(jTasFilt$attTypes == "factor")) {
                    message("Association Analysis : GLM")
                    message("(NOTE) Factors detected - running initial BLUE calculation...")
                    blueOut <- jRC$association(
                        rJava::.jnull(),
                        rJava::.jnull(),
                        jTasFilt$phenotype,
                        rJava::.jnull(),
                        as.integer(minClassSize),
                        biallelicOnly,
                        appendAddDom,
                        saveToFile,
                        outputFile,
                        maxP
                    )
                    blueOut <- blueOut$get("BLUE")
                    message("(NOTE) BLUEs calculated - using output to test markers...")
                    blueOut <- combineTasselGenotypePhenotype(
                        genotypeTable = tasObj@jGenotypeTable,
                        phenotype = blueOut
                    )
                    assocOut <- jRC$association(
                        rJava::.jnull(),
                        blueOut$genotypeTable(),
                        blueOut$phenotype(),
                        blueOut,
                        as.integer(minClassSize),
                        biallelicOnly,
                        appendAddDom,
                        saveToFile,
                        outputFile,
                        maxP
                    )
                    assocType <- "GLM"
                } else {
                    message("Association Analysis : GLM")
                    assocOut <- jRC$association(
                        rJava::.jnull(),
                        jTasFilt$genotypeTable,
                        jTasFilt$phenotype,
                        jTasFilt$genotypePhenotype,
                        as.integer(minClassSize),
                        biallelicOnly,
                        appendAddDom,
                        saveToFile,
                        outputFile,
                        maxP
                    )
                    assocType <- "GLM"
                }
            } else {
                stop("No TASSEL genotype table was found in TasselGenotypePhenotype object!")
            }
        } else {
            if (!rJava::is.jnull(tasObj@jGenotypeTable)) {
                if (any(apply(tmpDF, 2, function(x) any(is.na(x))))) {
                    stop("Missing phenotype data entries detected!")
                } else if (any(jTasFilt$attTypes == "factor")) {
                    message("Association Analysis : Fast Association")
                    message("(NOTE) Factors detected - running initial BLUE calculation...")
                    blueOut <- jRC$association(
                        rJava::.jnull(),
                        rJava::.jnull(),
                        jTasFilt$phenotype,
                        rJava::.jnull(),
                        as.integer(minClassSize),
                        biallelicOnly,
                        appendAddDom,
                        saveToFile,
                        outputFile,
                        maxP
                    )
                    blueOut <- blueOut$get("BLUE")
                    message("(NOTE) BLUEs calculated - using output to test markers...")
                    blueOut <- combineTasselGenotypePhenotype(
                        genotypeTable = tasObj@jGenotypeTable,
                        phenotype = blueOut
                    )
                    assocOut <- jRC$fastAssociation(
                        blueOut,
                        maxP,
                        maxThreads,
                        saveToFile,
                        outputFile
                    )
                    assocType <- "FastAssoc"
                } else {
                    message("Association Analysis : Fast Association")
                    assocOut <- jRC$fastAssociation(
                        jTasFilt$genotypePhenotype,
                        maxP,
                        maxThreads,
                        saveToFile,
                        outputFile
                    )
                    assocType <- "FastAssoc"
                }
            } else {
                stop("No TASSEL genotype table was found in TasselGenotypePhenotype object!")
            }
        }
    } else if (fitMarkers & !is.null(kinship) & !fastAssociation) {
        if (!rJava::is.jnull(tasObj@jGenotypeTable)) {
            if (any(jTasFilt$attTypes == "factor")) {
                message("Association Analysis : MLM")
                message("(NOTE) Factors detected - running initial BLUE calculation...")
                blueOut <- jRC$association(
                    rJava::.jnull(),
                    rJava::.jnull(),
                    jTasFilt$phenotype,
                    rJava::.jnull(),
                    as.integer(minClassSize),
                    biallelicOnly,
                    appendAddDom,
                    saveToFile,
                    outputFile,
                    maxP
                )
                blueOut <- blueOut$get("BLUE")
                message("(NOTE) BLUEs calculated - using output to test markers...")
                blueOut <- combineTasselGenotypePhenotype(
                    genotypeTable = tasObj@jGenotypeTable,
                    phenotype = blueOut
                )
                assocOut <- jRC$association(
                    kinship,
                    blueOut$genotypeTable(),
                    blueOut$phenotype(),
                    blueOut,
                    as.integer(minClassSize),
                    biallelicOnly,
                    appendAddDom,
                    saveToFile,
                    outputFile,
                    maxP
                )
                assocType <- "MLM"
            } else {
                message("Association Analysis : MLM")
                assocOut <- jRC$association(
                    kinship,
                    jTasFilt$genotypeTable,
                    jTasFilt$phenotype,
                    jTasFilt$genotypePhenotype,
                    as.integer(minClassSize),
                    biallelicOnly,
                    appendAddDom,
                    saveToFile,
                    outputFile,
                    maxP
                )
                assocType <- "MLM"
            }
        } else {
            stop("No TASSEL genotype table was found in TasselGenotypePhenotype object!")
        }
    } else {
        stop("Don't know how to analyze with given parameter inputs.")
    }

    return(tableReportList(assocOut))
}



## Phenotype filter - return modified TASSEL object - not exported (house keeping)
#' @importFrom rlang .data
tasPhenoFilter <- function(tasObj, filtObj) {

    # Get all TASSEL object trait metadata
    phenoAttDf <- extractPhenotypeAttDf(tasObj@jPhenotypeTable)

    # Get phenotype data frame
    phenoDF <- as.data.frame(tableReportToDF(tasObj@jPhenotypeTable))

    # Convert <data> and <covariates> to doubles (correct pass to TASSEL)
    doubCols <- as.character(
        phenoAttDf$traitName[which(phenoAttDf$traitType == "data" | phenoAttDf$traitType == "covariate")]
    )
    phenoDF[doubCols] <- sapply(phenoDF[doubCols], as.double)

    # Get taxa column
    taxaCol <- as.character(phenoAttDf$traitName[which(phenoAttDf$traitType == "taxa")])
    taxaNames <- as.vector(phenoDF[[taxaCol]])

    # Get non-taxa columns and reorder filtered columns (correct pass to TASSEL)
    origColNames <- colnames(phenoDF)
    filtObjRight <- c(taxaCol, filtObj)
    filtObjRight <- filtObjRight[match(origColNames, filtObjRight)]
    filtObjRight <- filtObjRight[!is.na(filtObjRight)]

    # Filter data frame columns based on association formula
    phenoDF <- phenoDF[, filtObjRight]
    phenoAttDf <- subset(phenoAttDf, subset = traitName %in% filtObjRight)

    # Get vector of non-taxa column names
    phenoColNames <- colnames(phenoDF)
    notTaxaCols <- phenoColNames[!(phenoColNames %in% taxaCol)]

    # Get attribute types
    attTypes <- as.vector(phenoAttDf$traitType[which(phenoAttDf$traitType != "taxa")])

    # Send filtered data frame to TASSEL methods
    jList <- rJava::new(rJava::J("java/util/ArrayList"))
    for (col_i in notTaxaCols) {
        jList$add(rJava::.jarray(phenoDF[[col_i]]))
    }

    # Create filtered TASSEL phenotype
    jc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    jc <- jc$createPhenotypeFromRDataFrameElements(
        taxaNames,
        rJava::.jarray(notTaxaCols),
        rJava::.jarray(attTypes),
        jList
    )

    # Return modified TASSEL objects
    if (rJava::is.jnull(tasObj@jGenotypeTable)) {
        return(
            list(
                attTypes = attTypes,
                phenoDf = phenoDF,
                finalVars = filtObj,
                notTaxaCols = notTaxaCols,
                genotypeTable = rJava::.jnull(),
                phenotype = jc,
                genotypePhenotype = rJava::.jnull()
            )
        )
    } else {
        jcComb <- combineTasselGenotypePhenotype(
            genotypeTable = tasObj@jGenotypeTable,
            phenotype = jc
        )
        return(
            list(
                attTypes = attTypes,
                phenoDf = phenoDF,
                finalVars = filtObj,
                notTaxaCols = notTaxaCols,
                genotypeTable = jcComb$genotypeTable(),
                phenotype = jcComb$phenotype(),
                genotypePhenotype = jcComb
            )
        )
    }
}


