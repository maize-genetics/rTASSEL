## ----
#' @title Plot SNP density across chromosomes
#'
#' @description Generates a heatmap-style visualization of SNP density
#'    across all chromosomes in a genotype table. Each chromosome is
#'    displayed on the Y-axis with genomic position on the X-axis. The
#'    number of SNPs within each window determines tile color intensity
#'    using the Viridis palette.
#'
#' @param tasObj An object of class \code{TasselGenotypePhenotype} or
#'    \code{TasselGenotype} that contains a genotype table.
#' @param windowSize Size of the genomic window (in base pairs) used to
#'    bin SNPs for density calculation. Defaults to \code{1e6} (1 Mbp).
#' @param colorOption Which viridis color palette to use? Options are:
#'    \code{"viridis"} (default), \code{"magma"}, \code{"inferno"},
#'    \code{"plasma"}, \code{"cividis"}, \code{"rocket"}, \code{"mako"},
#'    and \code{"turbo"}.
#' @param logNorm Should SNP counts be log\eqn{_{10}}{10}-transformed before
#'    mapping to fill color? Useful when a few windows have very high counts
#'    that compress the color scale. Defaults to \code{FALSE}.
#' @param interactive Do you want to produce an interactive visualization?
#'    Defaults to \code{FALSE}.
#'
#' @return Returns a \code{ggplot2} object or a \code{plotly} object if
#'    \code{interactive = TRUE}.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava J
#' @importFrom rlang .data
#'
#' @export
plotSnpDensity <- function(
    tasObj,
    windowSize = 1e6,
    colorOption = c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"),
    logNorm = FALSE,
    interactive = FALSE
) {
    isTGP <- is(tasObj, "TasselGenotypePhenotype")
    isTG  <- is(tasObj, "TasselGenotype")

    if (!isTGP && !isTG) {
        rlang::abort("tasObj must be of class \"TasselGenotypePhenotype\" or \"TasselGenotype\"")
    }

    if (isTGP && rJava::is.jnull(tasObj@jGenotypeTable)) {
        rlang::abort("tasObj does not contain a Genotype object")
    }

    if (!is.numeric(windowSize) || windowSize <= 0) {
        rlang::abort("`windowSize` must be a positive numeric value")
    }

    colorOption <- rlang::arg_match(colorOption)

    if (!is.logical(logNorm) || length(logNorm) != 1) {
        rlang::abort("`logNorm` must be a single logical value")
    }

    pSnpDensityParams <- list(
        "tasObj"      = tasObj,
        "windowSize"  = windowSize,
        "colorOption" = colorOption,
        "logNorm"     = logNorm
    )

    if (!interactive) {
        plotSnpDensityCore(pSnpDensityParams)
    } else {
        plotSnpDensityCoreInteractive(pSnpDensityParams)
    }
}


## ----
# @title Prepare SNP density data from a genotype table
# @param params A list of parameter variables
# @importFrom rJava J
primeSnpDensityData <- function(params) {
    tasObj     <- params$tasObj
    windowSize <- params$windowSize

    if (is(tasObj, "TasselGenotype")) {
        jtsPL <- tasObj@jRefObj$positions()
    } else {
        jtsPL <- getPositionList(tasObj)
    }
    posData <- rJava::J(
        "net/maizegenetics/plugindef/GenerateRCode"
    )$genotypeTableToPositionListOfArrays(jtsPL)

    chrVec <- posData$chromosomes
    posVec <- posData$startPos

    posByChr <- split(posVec, chrVec)

    dfList <- lapply(names(posByChr), function(chr) {
        bins    <- floor(posByChr[[chr]] / windowSize)
        allBins <- seq(0, max(bins))
        snpCount <- tabulate(match(bins, allBins), nbins = length(allBins))
        data.frame(
            snpCount    = snpCount,
            windowStart = allBins * windowSize,
            windowMid   = allBins * windowSize + (windowSize / 2),
            Chr         = chr,
            stringsAsFactors = FALSE
        )
    })

    counts <- do.call(rbind, dfList)

    chromNames <- unique(chrVec)
    chromOrder <- chromNames[order(nchar(chromNames), chromNames)]
    counts$Chr <- factor(counts$Chr, levels = chromOrder)

    return(counts)
}


## ----
# @title Core visual engine for SNP density plotting
# @param params A list of parameter variables
# @importFrom rlang .data
plotSnpDensityCore <- function(params) {
    densityDf   <- primeSnpDensityData(params)
    windowSize  <- params$windowSize
    colorOption <- params$colorOption
    logNorm     <- params$logNorm

    maxPos <- max(densityDf$windowMid)
    if (maxPos >= 1e6) {
        posScale  <- 1e6
        posUnit   <- "Mbp"
    } else if (maxPos >= 1e3) {
        posScale  <- 1e3
        posUnit   <- "kbp"
    } else {
        posScale  <- 1
        posUnit   <- "bp"
    }

    tileWidth   <- windowSize / posScale
    windowLabel <- if (windowSize >= 1e6) {
        paste0("Window size: ", windowSize / 1e6, " Mbp")
    } else if (windowSize >= 1e3) {
        paste0("Window size: ", windowSize / 1e3, " kbp")
    } else {
        paste0("Window size: ", windowSize, " bp")
    }

    densityDf$snpCount[densityDf$snpCount == 0] <- NA

    if (logNorm) {
        densityDf$snpCount <- log10(densityDf$snpCount)
    }

    fillLabel <- if (logNorm) expression(log[10](SNP~Count)) else "SNP Count"

    p <- ggplot2::ggplot(data = densityDf) +
        ggplot2::aes(
            x    = .data$windowMid / posScale,
            y    = .data$Chr,
            fill = .data$snpCount
        ) +
        ggplot2::geom_tile(
            width  = tileWidth,
            height = 0.8
        ) +
        ggplot2::scale_fill_viridis_c(
            name     = fillLabel,
            option   = colorOption,
            na.value = "grey95"
        ) +
        ggplot2::scale_y_discrete(limits = rev) +
        ggplot2::xlab(paste0("Position (", posUnit, ")")) +
        ggplot2::ylab("Chromosome") +
        ggplot2::labs(subtitle = windowLabel) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )

    return(p)
}


## ----
# @title Interactive SNP density plot
# @param params A list of parameter variables
plotSnpDensityCoreInteractive <- function(params) {
    p <- plotSnpDensityCore(params)
    return(plotly::ggplotly(p))
}


