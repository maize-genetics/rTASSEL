#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   GenomicPrediction.R
# Description:   General functions for running genomic prediction
# Author:        Brandon Monier
# Created:       2020-06-16 at 16:22:57
# Last Modified: 2020-06-16 at 16:23:52
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house all necessary
#    general functions and wrappers for TASSEL genomic prediction
#    capabilities.
#--------------------------------------------------------------------

#' @title R interface for TASSEL's genomic prediction capabilities
#'
#' @description This function acts as a front-end for TASSEL's genomic
#'   prediction functions. This analysis method uses gBLUP (genomic BLUP) to
#'   predict phenotypes from genotypes. It proceeds by fitting a mixed model
#'   that uses kinship to capture covariance between taxa. The mixed model can
#'   calculate BLUPs for taxa that do not have phenotypes based on the
#'   phenotypes of lines with relationship information.
#'
#'   A phenotype dataset and a kinship matrix must be supplied as input to the
#'   method by selecting both then choosing Analysis/Genomic Selection. In
#'   addition to trait values, the phenotype dataset may also contain factors
#'   or covariates which will be used as fixed effects in the model. All taxa
#'   in the phenotype dataset can only appear once. No repeated values are
#'   allowed for a single taxon. When the analysis is run, the user is
#'   presented with the choice to run k-fold cross-validation. If cross-
#'   validation is selected, then the number of folds and the number of
#'   iterations can be entered. For each iteration and each fold within an
#'   iteration, the correlation between the observed and predicted values will
#'   be reported. If cross-validation is not selected, then the original
#'   observations, predicted values and PEVs (prediction error variance) will
#'   be reported for all taxa in the dataset.
#'
#'   When k-fold cross-validation is performed, only taxa with phenotypes and
#'   rows in the kinship matrix are used. That set of taxa are divided into k
#'   subsets of equal size. Each subset in turn is used as the validation set.
#'   Phenotypes of the individuals in the validation are set to 0 then
#'   predicted using the remaining individuals as the training set. The
#'   correlation (r) of the observed values and predicted values is calculated
#'   for the validation set and reported. The mean and standard deviation of
#'   the mean of the r's are calculated for each trait and reported in the
#'   comments section of the "Accuracy" data set that is output by the
#'   analysis. In general, the results are not very sensitive to the choice of
#'   k. The number of iterations affects the standard error of the mean for the
#'   accuracy estimates. The defaults of k = 5 and iterations = 20 will be
#'   adequate for most users.
#'
#' @name genomicPrediction
#' @rdname genomicPrediction
#'
#' @param tasPhenoObj An object of class \code{TasselGenotypePenotype} that
#'   contains a phenotype object.
#' @param kinship A TASSEL kinship object.
#' @param doCV Do you want to perform k-fold cross-validation? Defaults to
#'   \code{FALSE}.
#' @param kFolds Number of folds to be entered.
#' @param nIter Number of iterations to be ran.
#'
#' @return Returns a \code{DataFrame}-based data frame
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava J
#' @export
genomicPrediction <- function(tasPhenoObj, kinship, doCV = FALSE, kFolds, nIter) {
    ## Check for correct rTASSEL class
    if (class(tasPhenoObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    ## Check to see if rTASSEL class object contains phenotype table
    jGenoTable <- getPhenotypeTable(tasPhenoObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL phenotype object not found")
    }

    ## Get phenotype pointer object
    tasPhenoObj <- tasPhenoObj@jPhenotypeTable

    ## Check to see if kinship parameter is of rJava and DistanceMatrix class
    if (class(kinship) == "rJava" && kinship != "net/maizegenetics/taxa/distance/DistanceMatrix") {
        stop("TASSEL kinship object is not of DistanceMatrix class")
    }

    ## Check to see if doCV parameter and subsequent parameters are right
    if (!doCV) {
        kFolds <- 0L
        nIter  <- 0L
    } else if (doCV) {
        kFolds <- as.integer(kFolds)
        nIter  <- as.integer(nIter)
    }

    ## Run genomic selection
    plugin <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    genSelRes <- plugin$genomicSelection(tasPhenoObj, kinship, doCV, kFolds, nIter)
    genSelRes <- tableReportToDF(genSelRes)

    ## Convert columns to appropriate data types
    if (!doCV) {
        colsNum <- c("Observed", "Predicted", "PEV")
        genSelRes[colsNum] <- sapply(genSelRes[colsNum], as.numeric)
    } else if (doCV) {
        colsNum <- c("Iteration", "Fold", "Accuracy")
        genSelRes[colsNum] <- sapply(genSelRes[colsNum], as.numeric)
    }

    ## Return object
    return(genSelRes)
}


