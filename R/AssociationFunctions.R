#--------------------------------------------------------------------
# Script Name:   AssociationFunctions.R
# Description:   General functions for running association analysis
# Author:        Brandon Monier
# Created:       2019-03-29 at 13:51:02
# Last Modified: 2019-04-01 at 17:27:39
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house all necessary
#    general functions and wrappers for TASSEL association analyses.
#--------------------------------------------------------------------

# Association analysis front-end
assocModelFitter <- function(tasObj,
                             formula,
                             fitMarkers = FALSE,
                             kinship = NULL) {

    # Logic - Check for TasselGenotypePhenotype class
    if (!is(tasObj, "TasselGenotypePhenotype")) {
        stop("tasObj is not of class \"TasselGenotypePhenotype\"")
    }

    # Extract formula response and prediction components
    formResp <- all.vars(formula[[2]])
    formPred <- all.vars(formula[[3]])

    # Get all TASSEL object trait metadata
    jtsPheno <- rTASSEL:::getPhenotypeTable(tasObj)
    phenoAttDf <- rTASSEL:::extractPhenotypeAttDf(jtsPheno)

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
        message("Running all non <data> traits...")
        finalResp <- formResp
        finalPred <- as.vector(tasPred[which(tasPred$traitName != "."), ]$traitName)
    } else {
        finalResp <- formResp
        finalPred <- formPred
    }

    # Logic - Handle association analyses
    if (!fitMarkers & is.null(kinship)) {
        if (!is.jnull(tasObj@jPhenotypeTable)) {
            message("Association Analysis : BLUEs")
        } else {
            stop("No TASSEL phenotype table was found in TasselGenotypePhenotype object!")
        }
    } else if (fitMarkers & is.null(kinship)) {
        if (!is.jnull(tasObj@jGenotypeTable)) {
            message("Association Analysis : GLM")
        } else {
            stop("No TASSEL genotype table was found in TasselGenotypePhenotype object!")
        }
    } else if (fitMarkers & !is.null(kinship)) {
        if (!is.jnull(tasObj@jGenotypeTable)) {
            message("Association Analysis : MLM")
        } else {
            stop("No TASSEL genotype table was found in TasselGenotypePhenotype object!")
        }
    } else {
        stop("Don't know how to analyze with given parameter inputs.")
    }

    # DEBUG
    # return(
    #     list(
    #         formResp = formResp,
    #         formPred = formPred,
    #         tasResp = tasResp,
    #         tasPred = tasPred
    #     )
    # )
    tasPhenoFilter(
        tasObj = tasObj,
        filtObj = c(finalResp, finalPred)
    )

}


# TASSEL Table to R data frame converter
tasTableConvert <- function(stringTab) {
    obj <- unlist(strsplit(stringTab, split = "\n"))
    obj <- strsplit(obj, split = "\t")
    obj <- t(simplify2array(obj))
    colnames(obj) <- as.character(unlist(obj[1, ]))
    obj <- obj[-1, ]
    tibble::as_tibble(obj)
}


# Phenotype filter
tasPhenoFilter <- function(tasObj, filtObj, attributeTypes = NULL) {

    # Get all TASSEL object trait metadata
    phenoAttDf <- rTASSEL:::extractPhenotypeAttDf(tasObj@jPhenotypeTable)

    # Get phenotype data frame
    phenoDf <- tasObj@jPhenotypeTable
    phenoDf <- tasTableConvert(phenoDf$toStringTabDelim())

    # Get taxa column
    taxaCol <- as.character(phenoAttDf$traitName[which(phenoAttDf$traitType == "taxa")])
    taxaNames <- as.vector(phenoDf[[taxaCol]])

    # Filter data frame columns based on association formula
    phenoDf <- phenoDf[, c(taxaCol, filtObj)]
    phenoAttDf <- subset(phenoAttDf, subset = traitName %in% c(taxaCol, filtObj))

    # Get vector of non-taxa column names
    phenoColNames <- colnames(phenoDf)
    notTaxaCols <- phenoColNames[!(phenoColNames %in% taxaCol)]

    # Get attribute types
    attTypes <- as.vector(phenoAttDf$traitType[which(phenoAttDf$traitType != "taxa")])

    # Send filtered data frame to TASSEL methods
    jList <- rJava::new(rJava::J("java/util/ArrayList"))
    for (col_i in notTaxaCols) {
        jList$add(rJava::.jarray(phenoDf[[col_i]]))
    }
    jc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    jc <- jc$createPhenotypeFromRDataFrameElements(
        taxaNames,
        notTaxaCols,
        attTypes,
        jList
    )
    rTASSEL:::.tasselObjectConstructor(jc)
    # return(tasTableConvert(jc$toStringTabDelim()))
}
