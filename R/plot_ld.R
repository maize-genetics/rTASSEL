## ----
haploviewColors <- function(dprime, lod) {
    colors <- rep("#FFFFFF", length(dprime))
    validD  <- !is.na(dprime)
    validL  <- !is.na(lod)
    dOne    <- validD & dprime >= (1 - 1e-6)
    lodHigh <- validL & lod >= 2

    colors[dOne & !lodHigh] <- "#bebfeb"
    colors[dOne & lodHigh]  <- "#db0003"

    pinkMask <- lodHigh & !dOne & validD
    d <- pmin(1, pmax(0, dprime[pinkMask]))
    colors[pinkMask] <- grDevices::rgb(1, 1 - d, 1 - d)

    colors
}


## ----
#' @title Linkage disequilibrium plot
#'
#' @description Calculates linkage disequilibrium (LD) and generates
#'   a static plot using \code{ggplot2} graphics.
#'
#' @name plotLD
#' @rdname plotLD
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}. That
#'   contains a genotype table.
#' @param ldType How do you want LD calculated? Currently, the available
#'   options are \code{"All"} and \code{"SlidingWindow"}. If
#'   \code{All} is selected, LD will be calculated for every
#'   combination of sites in the alignment (NOTE: this may produce a
#'   massive series of combinations; use only on heavily filtered
#'   genotype tables). If \code{SlidingWindow} is selected, LD will
#'   be calculated for sites within a window of sites surrounding the
#'   current site. Defaults to \code{"All"}.
#' @param windowSize What size do you want your LD analysis window? If you
#'   have chosen \code{SlidingWindow} for the \code{ldType} parameter, you
#'   will need to specify window size.
#' @param hetCalls How should heterozygous calls be handled? Current options
#'   are \code{"ignore"} (ignore altogether), \code{"missing"}
#'   (set to missing), and \code{"third"} (treat as third state).
#' @param plotVal What LD value do you want to plot? Options are:
#'   \itemize{
#'     \item \code{r2}: \eqn{r^{2}} (Default parameter)
#'     \item \code{DPrime}: \eqn{D'}
#'     \item \code{pDiseq}: \emph{p}-value
#'   }
#' @param colorScheme Color palette for the heatmap cells. \code{"haploview"}
#'   uses the classic Haploview scheme where cell color is determined by
#'   both LOD score and \eqn{D'}: LOD < 2 and \eqn{D'} < 1 gives white,
#'   LOD < 2 and \eqn{D'} = 1 gives blue, LOD \eqn{\ge} 2 and
#'   \eqn{D'} < 1 gives shades of pink/red scaled by \eqn{D'}, and
#'   LOD \eqn{\ge} 2 and \eqn{D'} = 1 gives bright red. All other
#'   options use continuous viridis-family palettes: \code{"viridis"}
#'   (default), \code{"magma"}, \code{"inferno"}, \code{"plasma"},
#'   \code{"cividis"}, \code{"rocket"}, \code{"mako"}, and
#'   \code{"turbo"}.
#' @param ldBlocks An optional \code{GRanges} object specifying genomic
#'   regions to highlight as LD blocks on the plot. Each range defines a
#'   chromosome and start/end positions; sites falling within each range
#'   are outlined with a triangular border. An optional \code{label}
#'   metadata column (in \code{mcols}) adds text annotations to each
#'   block. Defaults to \code{NULL}.
#' @param genomicTrack Logical. If \code{TRUE}, a horizontal genomic track
#'   is drawn above the LD plot showing the approximate physical positions
#'   of each site along the chromosome. Site IDs are placed at their
#'   physical positions on the track, and line segments connect each
#'   physical position to the corresponding evenly-spaced position above
#'   the LD triangle. When data span multiple chromosomes, positions are
#'   displayed cumulatively. Defaults to \code{FALSE}.
#' @param verbose Display messages? Defaults to \code{TRUE}.
#'
#' @details Linkage disequilibrium between any set of polymorphisms can be
#'   estimated by initially filtering a genotype dataset and then using
#'   this function. At this time, \eqn{D'}, \eqn{r^{2}} and P-values will be estimated. The
#'   current version calculates LD between haplotypes with known phase only
#'   (unphased diploid genotypes are not supported; see PowerMarker or
#'   Arlequin for genotype support).
#'   \itemize{
#'     \item \eqn{D'} is the standardized disequilibrium coefficient, a useful
#'     statistic for determining whether recombination or homoplasy has
#'     occurred between a pair of alleles.
#'     \item \eqn{r^{2}} represents the correlation between alleles at two loci, which
#'     is informative for evaluating the resolution of association approaches.
#'   }
#'   \eqn{D'} and \eqn{r^{2}} can be calculated when only two alleles are present. If more
#'   than two alleles, only the two most frequent alleles are used. P-values
#'   are determined by a two-sided Fisher's Exact test is calculated. Since LD
#'   is meaningless when scored with very small sample sizes, a minimum of 20
#'   taxa must be present to calculate LD and there must be 2 or more minor
#'   alleles.
#'
#' @seealso \code{\link{linkageDiseq}}
#'
#' @return Returns a \code{ggplot2} object.
#'
#' @importFrom GenomicRanges seqnames start end
#' @importFrom grDevices rgb
#' @importFrom rlang .data
#'
#' @export
plotLD <- function(
    tasObj,
    ldType      = c("All", "SlidingWindow"),
    windowSize  = NULL,
    hetCalls    = c("missing", "ignore", "third"),
    plotVal     = c("r2", "DPrime", "pDiseq"),
    colorScheme = c("viridis", "magma", "inferno", "plasma",
                    "cividis", "rocket", "mako", "turbo", "haploview"),
    ldBlocks     = NULL,
    genomicTrack = FALSE,
    verbose      = TRUE
) {
    angle <- 135

    plotVal     <- match.arg(plotVal)
    ldType      <- match.arg(ldType)
    colorScheme <- match.arg(colorScheme)
    useHaploview <- colorScheme == "haploview"

    ldDF <- linkageDiseq(
        tasObj     = tasObj,
        ldType     = ldType,
        windowSize = windowSize,
        hetCalls   = hetCalls,
        verbose    = verbose
    )

    if (plotVal == "r2") plotVal <- "R^2"
    ldDF$coord1 <- paste0(ldDF$Locus1, "_", ldDF$Position1)
    ldDF$coord2 <- paste0(ldDF$Locus2, "_", ldDF$Position2)

    # Sort sites by chromosome then position numerically
    sites <- unique(data.frame(
        coord = c(ldDF$coord1, ldDF$coord2),
        locus = c(ldDF$Locus1, ldDF$Locus2),
        pos   = as.numeric(c(ldDF$Position1, ldDF$Position2)),
        stringsAsFactors = FALSE
    ))
    locusAsNum <- suppressWarnings(as.numeric(sites$locus))
    if (all(!is.na(locusAsNum))) {
        sites <- sites[order(locusAsNum, sites$pos), ]
    } else {
        sites <- sites[order(sites$locus, sites$pos), ]
    }
    rownames(sites) <- NULL
    ids <- sites$coord

    # Reorder rows to match lower-triangle traversal in genomic order
    siteRank <- setNames(seq_len(nrow(sites)), sites$coord)
    r1 <- siteRank[ldDF$coord1]
    r2 <- siteRank[ldDF$coord2]
    swap <- r1 < r2
    if (any(swap)) {
        tmp <- ldDF$coord1[swap]
        ldDF$coord1[swap] <- ldDF$coord2[swap]
        ldDF$coord2[swap] <- tmp
    }
    ldDF <- ldDF[order(pmax(r1, r2), pmin(r1, r2)), ]

    if (useHaploview) {
        lod <- ldDF$N * ldDF[["R^2"]] / (2 * log(10))
        cellColors <- haploviewColors(ldDF$DPrime, lod)
        ldSub <- data.frame(
            coord1 = ldDF$coord1,
            coord2 = ldDF$coord2,
            val    = cellColors,
            stringsAsFactors = FALSE
        )
    } else {
        ldSub <- ldDF[, c("coord1", "coord2", plotVal)]
    }
    ldSub    <- as.data.frame(ldSub)
    ldSubRot <- ldCellRotater(ldSub, angle)

    nSites <- length(ids)
    idxLabels <- seq_len(nSites)

    idCoord <- list(
        "x" = seq(1.5, 1.5 + nSites - 1, 1),
        "y" = seq(0.5, 0.5 + nSites - 1, 1)
    )
    idCoord <- rotate(idCoord$x, idCoord$y, angle)

    scaleFactor   <- sqrt(nSites / 15)
    labelGap      <- 0.15 * scaleFactor
    borderLw      <- min(1, 10 / nSites)
    idxTextSize   <- max(2, min(5, 55 / nSites))
    snpTextSize   <- max(2, min(4, 60 / nSites))
    blockTextSize <- max(2.5, min(6, 65 / nSites))

    idxLabelHeight <- idxTextSize * 0.2 * scaleFactor
    idxLabelTop    <- max(idCoord$y) + labelGap + idxLabelHeight

    # Keep SNP-ID-to-track spacing stable across genomic windows.
    # Using observed max label length here causes padding to inflate as
    # coordinates get more digits in larger view spaces.
    refSnpNchar <- 10
    snpIdHeight <- refSnpNchar * snpTextSize * 0.15 * scaleFactor

    if (!is.null(ldBlocks)) {
        blockLabelSpace <- 0.3 * scaleFactor
        blockTop <- idxLabelTop + blockLabelSpace
        snpIdY <- blockTop + labelGap
    } else {
        snpIdY <- idxLabelTop + labelGap
    }

    upperMarginPadding <- snpIdY + snpIdHeight

    if (genomicTrack) {
        xRange <- range(idCoord$x)
        chroms <- unique(sites$locus)

        if (length(chroms) == 1) {
            posVals <- sites$pos
        } else {
            if (verbose) {
                message(
                    "Note: genomic track spans ", length(chroms),
                    " chromosomes; positions shown cumulatively."
                )
            }
            posVals   <- numeric(nrow(sites))
            cumOffset <- 0
            for (ch in chroms) {
                mask  <- sites$locus == ch
                chPos <- sites$pos[mask]
                posVals[mask] <- chPos - min(chPos) + cumOffset
                cumOffset <- max(posVals[mask]) + diff(range(chPos)) * 0.1 + 1
            }
        }

        posRange <- range(posVals)
        if (diff(posRange) > 0) {
            frac <- (posVals - posRange[1]) / diff(posRange)
        } else {
            frac <- rep(0.5, length(posVals))
        }
        physX <- xRange[2] - frac * diff(xRange)

        dataHeight         <- snpIdY - min(ldSubRot$y)
        snpIdTop           <- snpIdY + snpIdHeight * 0.55
        trackLineY         <- snpIdTop + dataHeight * 0.08
        trackTickSize      <- 0.1 * scaleFactor
        upperMarginPadding <- trackLineY + dataHeight * 0.02
    }

    if (!useHaploview) {
        legendLab <- switch(
            EXPR = plotVal,
            "R^2" = bquote(italic(r)^2),
            "DPrime" = bquote(italic(D)*"'"),
            "pDiseq" = bquote(italic(p)*'-value')
        )
    }

    cellPlot <- ggplot2::ggplot(data = ldSubRot) +
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$val, group = .data$group) +
        ggplot2::geom_polygon(color = "#D9D9D9", linewidth = borderLw) +
        ggplot2::annotate(
            geom = "text",
            x = idCoord$x,
            y = idCoord$y + labelGap,
            label = idxLabels,
            angle = 0,
            vjust = 0,
            size = idxTextSize
        )

    cellPlot <- cellPlot +
        ggplot2::annotate(
            geom = "text",
            x = idCoord$x,
            y = snpIdY,
            label = ids,
            angle = 90,
            hjust = 0,
            size = snpTextSize
        )

    if (genomicTrack) {
        cellPlot <- cellPlot +
            ggplot2::annotate(
                geom = "segment",
                x = idCoord$x, xend = physX,
                y = snpIdTop, yend = trackLineY,
                linewidth = 0.3, color = "gray50"
            ) +
            ggplot2::annotate(
                geom = "segment",
                x = min(idCoord$x) - 0.3, xend = max(idCoord$x) + 0.3,
                y = trackLineY, yend = trackLineY,
                linewidth = 0.6, color = "gray30"
            ) +
            ggplot2::annotate(
                geom = "segment",
                x = physX, xend = physX,
                y = trackLineY - trackTickSize,
                yend = trackLineY + trackTickSize,
                linewidth = 0.5, color = "gray30"
            )
    }

    cellPlot <- cellPlot +
        ggplot2::ylim(min(ldSubRot$y), upperMarginPadding) +
        ggplot2::scale_x_reverse() +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(
            legend.position = "bottom"
        )

    if (useHaploview) {
        cellPlot <- cellPlot + ggplot2::scale_fill_identity(na.value = "white")
    } else {
        cellPlot <- cellPlot +
            ggplot2::scale_fill_viridis_c(
                name     = legendLab,
                option   = colorScheme,
                na.value = "white"
            )
    }
    if (!is.null(ldBlocks)) {
        if (!inherits(ldBlocks, "GRanges")) {
            stop("'ldBlocks' must be a 'GRanges' object")
        }

        blockMcols <- S4Vectors::mcols(ldBlocks)
        hasLabels  <- "label" %in% colnames(blockMcols)
        blockLw    <- max(0.5, borderLw * 2)

        for (b in seq_len(length(ldBlocks))) {
            blockChr   <- as.character(GenomicRanges::seqnames(ldBlocks)[b])
            blockStart <- GenomicRanges::start(ldBlocks)[b]
            blockEnd   <- GenomicRanges::end(ldBlocks)[b]

            inBlock <- sites$locus == blockChr &
                sites$pos >= blockStart &
                sites$pos <= blockEnd

            if (!any(inBlock)) {
                warning("Block ", b, " contains no sites, skipping")
                next
            }

            blockIdx <- which(inBlock)
            s <- min(blockIdx)
            e <- max(blockIdx)

            blockCorners <- rotate(
                x = c(s, e + 1, s),
                y = c(s - 1, e, e),
                angle = angle
            )

            rectXmin <- min(blockCorners$x[1], blockCorners$x[2])
            rectXmax <- max(blockCorners$x[1], blockCorners$x[2])
            rectYmin <- max(blockCorners$y)

            combX <- c(rectXmin, rectXmin, blockCorners$x[3], rectXmax, rectXmax)
            combY <- c(blockTop, rectYmin, blockCorners$y[3], rectYmin, blockTop)

            cellPlot <- cellPlot +
                ggplot2::annotate(
                    geom = "polygon",
                    x = combX,
                    y = combY,
                    fill = NA,
                    color = "black",
                    linewidth = blockLw
                )

            if (hasLabels && !is.na(blockMcols$label[b])) {
                blockSizeBp <- blockEnd - blockStart
                if (blockSizeBp >= 1e6) {
                    sizeTag <- paste0("(", round(blockSizeBp / 1e6, 2), " Mbp)")
                } else {
                    sizeTag <- paste0("(", round(blockSizeBp / 1e3, 2), " kbp)")
                }
                blockLabel <- paste(blockMcols$label[b], sizeTag)

                labelX <- rectXmax - 0.1 * scaleFactor
                labelY <- idxLabelTop + blockLabelSpace / 3

                cellPlot <- cellPlot +
                    ggplot2::annotate(
                        geom = "text",
                        x = labelX,
                        y = labelY,
                        label = blockLabel,
                        hjust = 0,
                        size = blockTextSize,
                        fontface = "bold"
                    )
            }

            cellPlot <- cellPlot +
                ggplot2::annotate(
                    geom = "text",
                    x = idCoord$x[blockIdx],
                    y = snpIdY,
                    label = ids[blockIdx],
                    angle = 90,
                    hjust = 0,
                    size = snpTextSize,
                    fontface = "bold"
                )

            if (genomicTrack) {
                cellPlot <- cellPlot +
                    ggplot2::annotate(
                        geom = "segment",
                        x = idCoord$x[blockIdx], xend = physX[blockIdx],
                        y = snpIdTop, yend = trackLineY,
                        linewidth = 0.6, color = "black"
                    )
            }

        }
    }

    return(cellPlot)
}
