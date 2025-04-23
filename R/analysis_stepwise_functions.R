## ----
#' @title
#' Stepwise Model Fitter
#'
#' @description
#' Fits a model using stepwise selection criteria based on the provided
#' genotype-phenotype object.
#'
#' @param tasObj
#' A \code{TasselGenotypePhenotype} object. The input data for model 
#' fitting.
#' @param formula
#' A model formula, default is \code{. ~ .}.
#' @param modelType
#' Character. The model selection criteria used to determine which 
#' terms enter the model and how many. Must be one of \code{"pvalue"}, 
#' \code{"bic"}, \code{"mbic"}, or \code{"aic"}. Default is 
#' \code{"pvalue"}.
#' @param entryLimit
#' Numeric. The enter limit or maximum p-value for which a term can 
#' enter the model. Must be in range 0.0–1.0. Default is \code{0.01}.
#' @param exitLimit
#' Numeric. A term exits the model on a backward step if its p-value 
#' is greater than this value. Must be in range 0.0–1.0. Default is 
#' \code{0.01}.
#' @param maxNumberOfMarkers
#' Integer. The maximum number of markers that will be fit, if the 
#' enter limit is not reached first. Range 0–10000. Default is 
#' \code{20}.
#' @param nPermutations
#' Integer. Number of permutations for the model to determine an 
#' empirical alpha. Range 0–100000. Default is \code{0}.
#'
#' @examples
#' \dontrun{
#' results <- stepwiseModelFitter(tasObj)
#' }
#' 
#' @return
#' A list of association result tables including ANOVA reports and 
#' marker effect estimates, with and without confidence intervals.
#' 
#' @export
stepwiseModelFitter <- function(
    tasObj,
    formula = . ~ .,
    modelType = c("pvalue", "bic", "mbic", "aic"),
    entryLimit = 0.01,
    exitLimit = 0.01,
    maxNumberOfMarkers = 20,
    nPermutations = 0
) {
    # Logic - Check for TasselGenotypePhenotype class
    if (!is(tasObj, "TasselGenotypePhenotype")) {
        stop("tasObj is not of class \"TasselGenotypePhenotype\"")
    }

    # Logic - Check to see if TASSEL object has a phenotype table
    if (rJava::is.jnull(tasObj@jPhenotypeTable)) {
        stop("tasObj does not contain a Phenotype object")
    }

    # Logic - Check if formula is "all-by-all"
    if (formula[[2]] != "." || formula[[3]] != ".") {
        # Subset phenotype data
        rData        <- tableReportToDF(tasObj@jPhenotypeTable)
        attrData     <- makeAttributeData(tasObj@jPhenotypeTable, rData)
        traitsToKeep <- parseFormula(formula, attrData)
        subPh        <- selectTraitsFromJavaRef(tasObj@jPhenotypeTable, unlist(traitsToKeep))
        jSubPheno    <- javaRefObj(subPh)

        # Combine sub data with genotype
        jGtPh  <- combineTasselGenotypePhenotype(tasObj@jGenotypeTable, jSubPheno)
        tasObj <- .tasselObjectConstructor(jGtPh)
    }

    # Validate modelType
    stepPlugin <- rJava::.jnew(TASSEL_JVM$STEPWISE_PLUGIN)$modelType()
    modelType <- rlang::arg_match(modelType)

    modelType <- switch (modelType,
        "pvalue" = stepPlugin$pvalue,
        "bic"    = stepPlugin$bic,
        "mbic"   = stepPlugin$mbic,
        "aic"    = stepPlugin$aic
    )

    # Validate ranges
    if (!is.numeric(entryLimit) || entryLimit < 0 || entryLimit > 1) {
        rlang::abort("`entryLimit` must be a numeric value between 0 and 1.")
    }
    if (!is.numeric(exitLimit) || exitLimit < 0 || exitLimit > 1) {
        rlang::abort("`exitLimit` must be a numeric value between 0 and 1.")
    }
    if (!is.numeric(maxNumberOfMarkers) || maxNumberOfMarkers < 0 || maxNumberOfMarkers > 10000) {
        rlang::abort("`maxNumberOfMarkers` must be a numeric value between 0 and 10000.")
    }
    if (!is.numeric(nPermutations) || nPermutations < 0 || nPermutations > 100000) {
        rlang::abort("`nPermutations` must be a numeric value between 0 and 100000.")
    }

    modelFitter <- rJava::.jnew(
        TASSEL_JVM$STEPWISE_FITTER,
        tasObj@jTasselObj,
        "rt_stepwise"
    )

    # Set parameters
    modelFitter$setModelType(modelType)
    modelFitter$setEnterlimit(entryLimit)
    modelFitter$setExitlimit(exitLimit)
    modelFitter$setMaxNumberOfMarkers(as.integer(maxNumberOfMarkers))
    modelFitter$setNested(FALSE)
    modelFitter$setNumberOfPermutations(as.integer(nPermutations))

    # Run and report
    modelFitter$runAnalysis()

    stepAnovaReport   <- modelFitter$getAnovaReport()
    markerEstimates   <- modelFitter$getMarkerEffectReport()
    stepAnovaReportCi <- modelFitter$getAnovaReportWithCI()
    markerEstimatesCi <- modelFitter$getMarkerEffectReportWithCI()

    trl <- list(
        "ANOVA_report"        = tibble::as_tibble(tableReportToDF(stepAnovaReport)),
        "ANOVA_report_ci"     = tibble::as_tibble(tableReportToDF(stepAnovaReportCi)),
        "marker_estimates"    = tibble::as_tibble(tableReportToDF(markerEstimates)),
        "marker_estimates_ci" = tibble::as_tibble(tableReportToDF(markerEstimatesCi))
    )

    return(
        tableReportListToAssociationResults(
            trl = trl,
            aType = "Stepwise"
        )
    )
}


