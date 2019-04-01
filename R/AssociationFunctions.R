#--------------------------------------------------------------------
# Script Name:   AssociationFunctions.R
# Description:   General functions for running association analysis
# Author:        Brandon Monier
# Created:       2019-03-29 at 13:51:02
# Last Modified:
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

    # Get TASSEL object trait metadata
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
        stop("Variables in formula do not match traits in TASSEL object")
    } else if (any(!(formResp %in% tasResp$traitName))) {
        stop("Only <data> trait types can be implemented as response variables.")
    } else if (any(!(formPred %in% tasPred$traitName))) {
        stop("Only <factor> or <covariate> trait types can be implemented as predictor variables.")
    }

    # Logic - Handle "." variables
    if (formResp == "." & formPred == ".") {
        message("Running all traits...")
    } else if (formResp == "." & formPred != ".") {
        message("Running all <data> traits...")
    } else if (formResp != "." & formPred == ".") {
        message("Running all non <data> traits...")
    }

    # Logic - Handle association analyses
    if (!fitMarkers & is.null(kinship)) {
        message("Association Analysis : BLUEs")
    } else if (fitMarkers & is.null(kinship)) {
        if (!is.jnull(tasObj@jGenotypeTable)) {
            message("Association Analysis : GLM")
        } else {
            stop("No TASSEL genotype table was found in TasselGenotypePhenotype object!")
        }
    } else if (fitMarkers & !is.null(kinship)) {
        message("Associaiton Analysis : MLM")
    } else {
        stop("Don't know how to analyze with given parameter inputs.")
    }

    # DEBUG
    return(
        list(
            formResp = formResp,
            formPred = formPred,
            tasResp = tasResp,
            tasPred = tasPred
        )
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
