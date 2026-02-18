## ----
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
#' @importFrom rJava .jnull
#' @importFrom rJava new
#' @importFrom rJava J
#'
#' @export
# nocov start
treeJavaApp <- function(tasObj, clustMethod = c("Neighbor_Joining", "UPGMA")) {
    warnMsg <- paste0("The function 'ldJavaApp()' will be deprecated soon.")
    message(warnMsg)

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
# nocov end


