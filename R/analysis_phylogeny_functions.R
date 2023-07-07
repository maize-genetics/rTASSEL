#--------------------------------------------------------------------
# Script Name:   AnalysisPhylogenyFunctions.R
# Description:   General functions for running phylogeny analysis
# Author:        Brandon Monier
# Created:       2021-07-22 at 15:47:44
# Last Modified: 2021-07-22 at 16:09:29
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house all necessary
#    general functions and wrappers for TASSEL phylogeny analyses.
#--------------------------------------------------------------------

#' @title R interface for TASSEL's tree creation methods
#'
#' @description This function acts as a wrapper for TASSEL's
#'    \code{CreateTreePlugin}.
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#' @param clustMethod What clustering method should be used? Current options
#'    are \code{UGMA} and \code{Neighbor_Joining}. Defaults to
#'    \code{Neighbor_Joining}.
#'
#' @return Returns a \code{phylo} tree object. See the
#'    [ape]{https://cran.r-project.org/web/packages/ape/ape.pdf} package
#'    for further details.
#'
#' @importFrom ape read.tree
#' @importFrom rJava .jnull
#' @importFrom rJava new
#' @importFrom rJava J
#'
#' @export
createTree <- function(tasObj, clustMethod = c("Neighbor_Joining", "UPGMA")) {
    if (!is(tasObj, "TasselGenotypePhenotype")) {
        stop("tasObj is not of class \"TasselGenotypePhenotype\"")
    }

    clustMethod <- match.arg(clustMethod)

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

    ape::read.tree(text = myTree$toString(), tree.names = clustMethod)
}


