## ----
#' @title Generate PCA plot
#'
#' @description
#' This function will generate a general 2D scatterplot for a given set
#' of principal components.
#'
#' @param pcaObj A \code{PCAResults} object containing eigenvalue data (
#'    e.g \code{PC_Datum}).
#' @param x Principal component number to plot on x-axis.
#' @param y Principal component number to plot on y-axis.
#' @param cluster Plot points by hierachical clustering groups? Defaults to
#'    \code{FALSE}.
#' @param nClust Number of clusters to report.
#' @param metadata A \code{data.frame} object of additional categorical
#'    information for each sample/taxa. Defaults to \code{NULL}.
#' @param mCol What metadata column do you want to plot? Needs metadata
#'    object to work.
#' @param interactive Do you want to produce an interactive visualization?
#'    Defaults to \code{FALSE}.
#' @param pointColor color of non-categorical PCA data points.
#'
#' @export
plotPCA <- function(
    pcaObj,
    x = 1,
    y = 2,
    cluster = FALSE,
    nClust = 2,
    metadata = NULL,
    mCol = NULL,
    interactive = FALSE,
    pointColor = "black"
) {
    if (!is(pcaObj, "PCAResults")) {
        stop(
            "The object '", deparse(substitute(assocRes)),
            "' is not an 'AssociationResults' object"
        )
    }

    if (!"PC_Datum" %in% reportNames(pcaObj)) {
        stop("'PC_Datum' table report not found")
    }

    pPCACoreParams <- list(
        "pcaObj"     = pcaObj,
        "x"          = x,
        "y"          = y,
        "cluster"    = cluster,
        "nClust"     = nClust,
        "pointColor" = pointColor,
        "metaData"   = metadata,
        "mCol"       = mCol
    )

    if (!interactive) {
        plotPCACore(pPCACoreParams)
    } else {
        plotPCACoreInteractive(pPCACoreParams)
    }
}


## ----
# @title Core visual engine for PCA plotting
# @param params A list of parameter variables
# @importFrom rlang .data
plotPCACore <- function(params) {
    ## Parse parameters
    pcaObj     <- params$pcaObj
    x          <- params$x
    y          <- params$y
    cluster    <- params$cluster
    nClust     <- params$nClust
    pointColor <- params$pointColor
    metaData   <- params$metaData
    mCol       <- params$mCol

    pcDf <- tableReport(pcaObj, "PC_Datum")

    if (cluster) {
        pcDf <- hClustPCAResults(pcDf, nClust)
    }

    if (!is.null(metaData)) {
        if (!"Taxa" %in% colnames(metaData)) {
            stop("Metadata is missing 'Taxa' column")
        }
        if (!mCol %in% colnames(metaData)) {
            stop("Metadata column missing from metadata")
        }
        if (!any(pcDf$Taxa %in% metaData$Taxa)) {
            stop("No samples match what is found in given metadata column")
        }
        pcDf <- merge(pcDf, metaData, by = "Taxa")
    }

    if (grepl("PC", x)) x <- suppressWarnings(as.numeric(gsub("PC", "", x)))
    if (grepl("PC", y)) y <- suppressWarnings(as.numeric(gsub("PC", "", y)))

    if (!paste0("PC", x) %in% colnames(pcDf) || is.na(x)) {
        stop("Principal component value for 'x' is not found in data")
    }
    if (!paste0("PC", y) %in% colnames(pcDf) || is.na(y)) {
        stop("Principal component value for 'y' is not found in data")
    }

    ## Filter based on PCs
    pcIdMatches <- colnames(pcDf)[colnames(pcDf) %in% paste0("PC", c(x, y))]
    if (cluster && is.null(metaData)) {
        pcDfSub <- pcDf[, c("Taxa", pcIdMatches, "Cluster")]
    } else if (!is.null(metaData) && !cluster) {
        pcDfSub <- pcDf[, c("Taxa", pcIdMatches, mCol)]
    } else if (!is.null(metaData) && cluster) {
        pcDfSub <- pcDf[, c("Taxa", pcIdMatches, mCol, "Cluster")]
    } else {
        pcDfSub <- pcDf[, c("Taxa", pcIdMatches)]
    }

    ## Get proportion of total variance
    if (!"Eigenvalues_Datum" %in% reportNames(pcaObj)) {
        message("Eigenvalues are not reported and will not be shown in axis text")
        xAxisText <- paste0("PC", x)
        yAxisText <- paste0("PC", y)
    } else {
        eigVal <- tableReport(pcaObj, "Eigenvalues_Datum")
        xProp <- eigVal[as.numeric(eigVal$PC) + 1 == x, "proportion_of_total"]
        yProp <- eigVal[as.numeric(eigVal$PC) + 1 == y, "proportion_of_total"]
        xAxisText <- paste0("PC", x, " (", round(xProp * 100, 2), "%)")
        yAxisText <- paste0("PC", y, " (", round(yProp * 100, 2), "%)")
    }

    ## Clustering aesthetics
    if (cluster && is.null(metaData)) {
        pcaAesMain <- ggplot2::aes(
            x     = .data[[paste0("PC", x)]],
            y     = .data[[paste0("PC", y)]],
            color = .data[["Cluster"]]
        )
        pcaGeomPoint <- ggplot2::geom_point()
    } else if (!is.null(metaData) && !cluster) {
        pcaAesMain <- ggplot2::aes(
            x     = .data[[paste0("PC", x)]],
            y     = .data[[paste0("PC", y)]],
            color = .data[[mCol]]
        )
        pcaGeomPoint <- ggplot2::geom_point()
    } else if (!is.null(metaData) && cluster) {
        pcaAesMain <- ggplot2::aes(
            x     = .data[[paste0("PC", x)]],
            y     = .data[[paste0("PC", y)]],
            color = .data[[mCol]],
            shape = .data[["Cluster"]]
        )
        pcaGeomPoint <- ggplot2::geom_point()

    } else {
        pcaAesMain <- ggplot2::aes(
            x = .data[[paste0("PC", x)]],
            y = .data[[paste0("PC", y)]]
        )
        pcaGeomPoint <- ggplot2::geom_point(color = pointColor)
    }

    p <- ggplot2::ggplot(data = pcDfSub) +
        pcaAesMain +
        pcaGeomPoint +
        ggplot2::xlab(xAxisText) +
        ggplot2::ylab(yAxisText) +
        ggplot2::theme_bw()

    return(p)
}


## ----
# @title Modify fancy quotes for plotly
# @param params A list of parameter variables
plotPCACoreInteractive <- function(params) {
    p <- plotPCACore(params)
    return(plotly::ggplotly(p))
}


## ----
hClustPCAResults <- function(p, nClust) {
    distMat     <- stats::dist(p[, 2:ncol(p)], method = "euclidean")
    hClustAvg   <- stats::hclust(distMat, method = "average")
    clustAssign <- stats::cutree(hClustAvg, k = nClust)

    p$Cluster <- as.character(clustAssign)

    return(p)
}


