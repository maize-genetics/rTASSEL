## ----
#' @title Generate scree plots
#'
#' @description
#' This function will generate a line plot of eigenvalues of PCs. These
#' plots can be used to determine the number of factors to retain for
#' downstream analyses (e.g GWAS).
#'
#' @param pcaObj A \code{PCAResults} object containing eigenvalue data (
#'    e.g \code{Eigenvalues_Datum}).
#' @param nComp Total number of principal components to plot. Defaults to
#'    \code{10}.
#' @param lineColor Color of scree line.
#' @param pointColor Color of scree points.
#' @param interactive Do you want to produce an interactive visualization?
#'    Defaults to \code{FALSE}.
#'
#' @export
plotScree <- function(
    pcaObj,
    nComp = 10,
    interactive = FALSE,
    lineColor = "grey",
    pointColor = "black"
) {
    if (!is(pcaObj, "PCAResults")) {
        stop(
            "The object '", deparse(substitute(assocRes)),
            "' is not an 'AssociationResults' object"
        )
    }

    if (!"Eigenvalues_Datum" %in% reportNames(pcaObj)) {
        stop("'Eigenvalues_Datum' table report not found")
    }

    pScreeCoreParams <- list(
        "pcaObj"      = pcaObj,
        "nComp"       = nComp,
        "lineColor"   = lineColor,
        "pointColor"  = pointColor
    )

    if (!interactive) {
        plotScreeCore(pScreeCoreParams)
    } else {
        plotScreeCoreInteractive(pScreeCoreParams)
    }
}


## ----
# @title Core visual engine for scree plotting
# @param params A list of parameter variables
#' @importFrom rlang .data
#' @importFrom utils head
plotScreeCore <- function(params) {
    ## Parse parameters
    pcaObj     <- params$pcaObj
    nComp      <- params$nComp
    lineColor  <- params$lineColor
    pointColor <- params$pointColor

    eigenValueDf <- tableReport(pcaObj, "Eigenvalues_Datum")

    if (nComp > nrow(eigenValueDf)) {
        stop("Number of principal components exceed what is reported")
    }

    eigenValueDfSub <- utils::head(eigenValueDf, n = nComp)
    eigenValueDfSub$PC <- as.numeric(eigenValueDfSub$PC) + 1

    p <- ggplot2::ggplot(data = eigenValueDfSub) +
        ggplot2::aes(x = .data$PC, y = .data$proportion_of_total) +
        ggplot2::geom_line(color = "grey") +
        ggplot2::geom_point(color = "black") +
        ggplot2::xlab("Principal components") +
        ggplot2::ylab("Proportion of total variance") +
        ggplot2::theme_bw()

    return(p)
}


## ----
# @title Modify fancy quotes for plotly
# @param params A list of parameter variables
plotScreeCoreInteractive <- function(params) {
    p <- plotScreeCore(params)
    return(plotly::ggplotly(p))
}


