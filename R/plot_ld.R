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
computeSizingParams <- function(nSites) {
    scaleFactor <- sqrt(nSites / 15)
    snpTextSize <- max(2, min(4, 60 / nSites))
    list(
        scaleFactor   = scaleFactor,
        labelGap      = 0.15 * scaleFactor,
        borderLw      = min(1, 10 / nSites),
        idxTextSize   = max(2, min(5, 55 / nSites)),
        snpTextSize   = snpTextSize,
        blockTextSize = max(2.5, min(6, 65 / nSites)),
        snpIdHeight   = 10 * snpTextSize * 0.15 * scaleFactor
    )
}


## ----
computePhysicalPositions <- function(sites, verbose) {
    chroms <- unique(sites$locus)

    if (length(chroms) == 1) return(sites$pos)

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
    posVals
}


## ----
resolveBlocks <- function(ldBlocks, sites) {
    if (inherits(ldBlocks, "GRanges")) {
        blockMcols <- S4Vectors::mcols(ldBlocks)
        hasLabels  <- "label" %in% colnames(blockMcols)
        n <- length(ldBlocks)
        data.frame(
            chr       = as.character(GenomicRanges::seqnames(ldBlocks)),
            start     = GenomicRanges::start(ldBlocks),
            end       = GenomicRanges::end(ldBlocks),
            label     = if (hasLabels) as.character(blockMcols$label) else NA_character_,
            color     = rep("black", n),
            linewidth = rep(NA_real_, n),
            showSpan  = rep(TRUE, n),
            stringsAsFactors = FALSE
        )
    } else {
        if (methods::is(ldBlocks, "LDRegion")) {
            ldBlocks <- list(ldBlocks)
        }
        allLDRegion <- vapply(ldBlocks, methods::is, logical(1), class2 = "LDRegion")
        if (!all(allLDRegion)) {
            stop(
                "'ldBlocks' must be a GRanges object, an LDRegion object, ",
                "or a list of LDRegion objects"
            )
        }
        chroms <- unique(sites$locus)
        if (length(chroms) != 1) {
            stop(
                "LDRegion blocks require single-chromosome data; ",
                "use a GRanges object for multi-chromosome block specification"
            )
        }
        data.frame(
            chr       = rep(chroms, length(ldBlocks)),
            start     = vapply(ldBlocks, slot, numeric(1), "start"),
            end       = vapply(ldBlocks, slot, numeric(1), "end"),
            label     = vapply(ldBlocks, slot, character(1), "label"),
            color     = vapply(ldBlocks, slot, character(1), "color"),
            linewidth = vapply(ldBlocks, slot, numeric(1), "linewidth"),
            showSpan  = vapply(ldBlocks, slot, logical(1), "showSpan"),
            stringsAsFactors = FALSE
        )
    }
}


## ----
#' @title Linkage disequilibrium plot
#'
#' @description Generates a static LD heatmap from a pre-computed
#'   \code{\linkS4class{LDResults}} object using \code{ggplot2} graphics.
#'
#' @name plotLD
#' @rdname plotLD
#'
#' @param ldObj An object of class \code{\linkS4class{LDResults}} as
#'   returned by \code{\link{linkageDiseq}}.
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
#' @param ldBlocks Optional specification of genomic regions to highlight
#'   as LD blocks on the plot. Accepted formats:
#'   \itemize{
#'     \item A \code{\linkS4class{LDRegion}} object (single block) or a
#'       \code{list} of \code{LDRegion} objects (multiple blocks). Each
#'       region carries \code{start}, \code{end}, an optional \code{label},
#'       and an outline \code{color}. The chromosome is inferred from the
#'       data (single-chromosome only).
#'     \item A \code{GRanges} object where each range defines a chromosome
#'       and start/end positions. An optional \code{label} metadata column
#'       (in \code{mcols}) adds text annotations. Block outlines default
#'       to black.
#'   }
#'   Sites falling within each region are outlined with a triangular
#'   border. Defaults to \code{NULL}.
#' @param genomicTrack Logical. If \code{TRUE}, a horizontal genomic track
#'   is drawn above the LD plot showing the approximate physical positions
#'   of each site along the chromosome. Site IDs are placed at their
#'   physical positions on the track, and line segments connect each
#'   physical position to the corresponding evenly-spaced position above
#'   the LD triangle. When data span multiple chromosomes, positions are
#'   displayed cumulatively. Defaults to \code{FALSE}.
#' @param showIndex Logical. If \code{TRUE} (default), numeric index labels
#'   (1, 2, 3, \ldots) are drawn along the diagonal of the LD triangle.
#'   Set to \code{FALSE} to hide them.
#' @param verbose Display messages? Defaults to \code{TRUE}.
#'
#' @details This function visualises a pre-computed LD result set.
#'   Use \code{\link{linkageDiseq}} to calculate LD first, then pass the
#'   resulting \code{\linkS4class{LDResults}} object here.
#'
#' @seealso \code{\link{linkageDiseq}}, \code{\linkS4class{LDResults}}
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom GenomicRanges seqnames start end
#' @importFrom grDevices rgb
#' @importFrom methods is slot
#' @importFrom rlang .data
#'
#' @export
plotLD <- function(
    ldObj,
    plotVal     = c("r2", "DPrime", "pDiseq"),
    colorScheme = c("viridis", "magma", "inferno", "plasma",
                    "cividis", "rocket", "mako", "turbo", "haploview"),
    ldBlocks     = NULL,
    genomicTrack = FALSE,
    showIndex    = TRUE,
    verbose      = TRUE
) {
    if (!methods::is(ldObj, "LDResults")) {
        stop("'ldObj' must be an object of class \"LDResults\"")
    }

    angle <- 135

    plotVal     <- match.arg(plotVal)
    colorScheme <- match.arg(colorScheme)
    useHaploview <- colorScheme == "haploview"

    if (!is.null(ldBlocks)) {
        validType <- inherits(ldBlocks, "GRanges") ||
            methods::is(ldBlocks, "LDRegion") ||
            (is.list(ldBlocks) && all(vapply(
                ldBlocks, methods::is, logical(1), class2 = "LDRegion"
            )))
        if (!validType) {
            stop(
                "'ldBlocks' must be a GRanges object, an LDRegion object, ",
                "or a list of LDRegion objects"
            )
        }
    }

    ldDF <- as.data.frame(ldObj@results)

    if (plotVal == "r2") plotVal <- "R^2"
    ldDF$coord1 <- paste0(ldDF$Locus1, "_", ldDF$Position1)
    ldDF$coord2 <- paste0(ldDF$Locus2, "_", ldDF$Position2)

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

    idCoord <- list(
        "x" = seq(1.5, 1.5 + nSites - 1, 1),
        "y" = seq(0.5, 0.5 + nSites - 1, 1)
    )
    idCoord <- rotate(idCoord$x, idCoord$y, angle)

    sizing      <- computeSizingParams(nSites)
    labelGap    <- sizing$labelGap
    borderLw    <- sizing$borderLw
    idxTextSize <- sizing$idxTextSize
    snpTextSize <- sizing$snpTextSize
    scaleFactor <- sizing$scaleFactor
    snpIdHeight <- sizing$snpIdHeight

    # -- Heatmap + index labels ----------------------------------------
    p <- ggplot2::ggplot(data = ldSubRot) +
        ggplot2::aes(
            x = .data$x, y = .data$y,
            fill = .data$val, group = .data$group
        ) +
        ggplot2::geom_polygon(color = "#D9D9D9", linewidth = borderLw)

    if (showIndex) {
        p <- p + ggplot2::annotate(
            geom  = "text",
            x     = idCoord$x,
            y     = idCoord$y + labelGap,
            label = seq_len(nSites),
            angle = 0,
            vjust = 0,
            size  = idxTextSize
        )
    }

    idxLabelHeight <- if (showIndex) idxTextSize * 0.2 * scaleFactor else 0
    idxLabelTop    <- max(idCoord$y) + labelGap + idxLabelHeight

    hasBlocks <- !is.null(ldBlocks)
    if (hasBlocks) {
        blockLabelSpace <- 0.3 * scaleFactor
        blockTop <- idxLabelTop + blockLabelSpace
    }

    # -- Fill scale ----------------------------------------------------
    if (!useHaploview) {
        legendLab <- switch(
            EXPR = plotVal,
            "R^2"    = bquote(italic(r)^2),
            "DPrime" = bquote(italic(D) * "'"),
            "pDiseq" = bquote(italic(p) * "-value")
        )
        p <- p + ggplot2::scale_fill_viridis_c(
            name     = legendLab,
            option   = colorScheme,
            na.value = "white"
        )
    } else {
        p <- p + ggplot2::scale_fill_identity(na.value = "white")
    }

    # -- LD block outlines and labels ----------------------------------
    blockIndices <- list()
    if (hasBlocks) {
        resolvedDf     <- resolveBlocks(ldBlocks, sites)
        defaultBlockLw <- max(0.5, borderLw * 2)
        blockTextSize  <- sizing$blockTextSize

        for (b in seq_len(nrow(resolvedDf))) {
            blockChr   <- resolvedDf$chr[b]
            blockStart <- resolvedDf$start[b]
            blockEnd   <- resolvedDf$end[b]
            blockColor <- resolvedDf$color[b]
            blockLw    <- if (is.na(resolvedDf$linewidth[b])) defaultBlockLw else resolvedDf$linewidth[b]

            inBlock <- sites$locus == blockChr &
                sites$pos >= blockStart &
                sites$pos <= blockEnd

            if (!any(inBlock)) {
                warning("Block ", b, " contains no sites, skipping")
                next
            }

            blockIdx <- which(inBlock)
            blockIndices[[length(blockIndices) + 1]] <- blockIdx
            s <- min(blockIdx)
            e <- max(blockIdx)

            blockCorners <- rotate(
                x     = c(s, e + 1, s),
                y     = c(s - 1, e, e),
                angle = angle
            )

            rectXmin <- min(blockCorners$x[1], blockCorners$x[2])
            rectXmax <- max(blockCorners$x[1], blockCorners$x[2])
            rectYmin <- max(blockCorners$y)

            combX <- c(rectXmin, rectXmin, blockCorners$x[3], rectXmax, rectXmax)
            combY <- c(blockTop, rectYmin, blockCorners$y[3], rectYmin, blockTop)

            p <- p + ggplot2::annotate(
                geom      = "polygon",
                x         = combX,
                y         = combY,
                fill      = NA,
                color     = blockColor,
                linewidth = blockLw
            )

            hasLabel   <- !is.na(resolvedDf$label[b])
            wantsSpan  <- resolvedDf$showSpan[b]
            if (hasLabel || wantsSpan) {
                blockSizeBp <- blockEnd - blockStart
                sizeTag <- if (blockSizeBp >= 1e6) {
                    paste0("(", round(blockSizeBp / 1e6, 2), " Mbp)")
                } else {
                    paste0("(", round(blockSizeBp / 1e3, 2), " kbp)")
                }

                blockLabel <- if (hasLabel && wantsSpan) {
                    paste(resolvedDf$label[b], sizeTag)
                } else if (hasLabel) {
                    resolvedDf$label[b]
                } else {
                    sizeTag
                }

                p <- p + ggplot2::annotate(
                    geom     = "text",
                    x        = rectXmax - 0.1 * scaleFactor,
                    y        = idxLabelTop + blockLabelSpace / 3,
                    label    = blockLabel,
                    hjust    = 0,
                    size     = blockTextSize,
                    fontface = "bold",
                    color    = blockColor
                )
            }
        }
    }

    # -- SNP ID labels -------------------------------------------------
    snpIdY <- if (hasBlocks) blockTop + labelGap else idxLabelTop + labelGap

    p <- p + ggplot2::annotate(
        geom  = "text",
        x     = idCoord$x,
        y     = snpIdY,
        label = ids,
        angle = 90,
        hjust = 0,
        size  = snpTextSize
    )

    if (length(blockIndices) > 0) {
        allBlockIdx <- unique(unlist(blockIndices))

        p <- p + ggplot2::annotate(
            geom     = "text",
            x        = idCoord$x[allBlockIdx],
            y        = snpIdY,
            label    = ids[allBlockIdx],
            angle    = 90,
            hjust    = 0,
            size     = snpTextSize,
            fontface = "bold"
        )
    }

    topY <- snpIdY + snpIdHeight

    # -- Genomic track -------------------------------------------------
    if (genomicTrack) {
        posVals  <- computePhysicalPositions(sites, verbose)
        xRange   <- range(idCoord$x)
        posRange <- range(posVals)
        frac <- if (diff(posRange) > 0) {
            (posVals - posRange[1]) / diff(posRange)
        } else {
            rep(0.5, length(posVals))
        }
        physX <- xRange[2] - frac * diff(xRange)

        trackGap      <- 0.2 * scaleFactor
        trackLineY    <- (topY + trackGap) * 1.5
        trackTickSize <- 0.2 * scaleFactor

        p <- p +
            ggplot2::annotate(
                geom      = "segment",
                x         = idCoord$x,
                xend      = physX,
                y         = topY,
                yend      = trackLineY,
                linewidth = 0.3,
                color     = "gray50"
            ) +
            ggplot2::annotate(
                geom      = "segment",
                x         = min(idCoord$x) - 0.3,
                xend      = max(idCoord$x) + 0.3,
                y         = trackLineY,
                yend      = trackLineY,
                linewidth = 0.6,
                color     = "gray30"
            ) +
            ggplot2::annotate(
                geom      = "segment",
                x         = physX,
                xend      = physX,
                y         = trackLineY - trackTickSize,
                yend      = trackLineY + trackTickSize,
                linewidth = 0.5,
                color     = "gray30"
            )

        if (length(blockIndices) > 0) {
            p <- p +
                ggplot2::annotate(
                    geom      = "segment",
                    x         = idCoord$x[allBlockIdx],
                    xend      = physX[allBlockIdx],
                    y         = topY,
                    yend      = trackLineY,
                    linewidth = 0.6,
                    color     = "black"
                )
        }

        chroms <- unique(sites$locus)
        formatRange <- function(lo, hi) {
            span <- hi - lo
            if (span >= 1e6 || hi >= 1e6) {
                paste0(round(lo / 1e6, 2), " \u2013 ", round(hi / 1e6, 2), " Mbp")
            } else {
                paste0(round(lo / 1e3, 2), " \u2013 ", round(hi / 1e3, 2), " kbp")
            }
        }
        if (length(chroms) == 1) {
            rawRange   <- range(sites$pos)
            trackLabel <- paste0("Chr ", chroms, ": ", formatRange(rawRange[1], rawRange[2]))
        } else {
            parts <- vapply(chroms, function(ch) {
                chPos <- sites$pos[sites$locus == ch]
                paste0("Chr ", ch, ": ", formatRange(min(chPos), max(chPos)))
            }, character(1))
            trackLabel <- paste(parts, collapse = " | ")
        }

        trackLabelSize <- max(2.5, min(4, 50 / nSites))
        trackLabelY    <- trackLineY + trackTickSize + 0.4 * scaleFactor

        p <- p + ggplot2::annotate(
            geom     = "text",
            x        = max(xRange) + 0.3,
            y        = trackLabelY,
            label    = trackLabel,
            hjust    = 0,
            size     = trackLabelSize,
            color    = "gray20",
            fontface = "bold"
        )

        topY <- trackLabelY + trackLabelSize * 0.15 * scaleFactor
    }

    # -- Coord, theme, and final assembly ------------------------------
    p +
        ggplot2::scale_x_reverse() +
        ggplot2::ylim(min(ldSubRot$y), topY) +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "bottom")
}
