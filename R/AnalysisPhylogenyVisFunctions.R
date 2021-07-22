#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   AnalysisPhylogenyVisFunctions.R
# Description:   General functions for visualizing phylogeny
# Author:        Brandon Monier
# Created:       2021-07-22 at 16:50:16
# Last Modified: 2021-07-22 at 17:08:53
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house all necessary
#    general functions and wrappers for visualizing TASSEL phylogeny
#--------------------------------------------------------------------

#' @title R interface for Archaeopteryx interactive tree viewer
#'
#' @description This function acts as a wrapper for TASSEL's
#'    interface to the Archaeopteryx Java tree Viewer.
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#' @param clustMethod What clustering method should be used? Current options
#'    are \code{UGMA} and \code{Neighbor_Joining}. Defaults to
#'    \code{Neighbor_Joining}.
#'
#' @return Returns a Java-based visualization application.
#'
#' @importFrom ape read.tree
#' @importFrom rJava .jnull
#' @importFrom rJava new
#' @importFrom rJava J
#'
#' @export
treeJavaApp <- function(tasObj, clustMethod = c("Neighbor_Joining", "UPGMA")) {
    if (!is(tasObj, "TasselGenotypePhenotype")) {
        stop("tasObj is not of class \"TasselGenotypePhenotype\"")
    }

    clustMethod <- match.arg(clustMethod)

    # Get TASSEL tree object
    plugin <- rJava::new(
        rJava::J("net/maizegenetics/analysis/tree/CreateTreePlugin"),
        rJava::.jnull("java/awt/Frame"),
        FALSE
    )

    input <- rJava::J("net/maizegenetics/plugindef/DataSet")
    input <- input$getDataSet(getGenotypeTable(tasObj))

    plugin$setParameter("clusteringMethod", clustMethod)
    plugin$setParameter("saveDistanceMatrix", "false")
    myTree <- plugin$runPlugin(input)


    # Call Archaeopteryx
    archPlugin <- rJava::new(
        rJava::J("net/maizegenetics/analysis/tree/ArchaeopteryxPlugin"),
        rJava::.jnew("java/awt/Frame"),
        TRUE
    )
    treeInput <- rJava::J("net/maizegenetics/plugindef/DataSet")
    treeInput <- input$getDataSet(myTree)
    archPlugin$performFunction(treeInput)
}


