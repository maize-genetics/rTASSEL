## ----
# @title Build zoom connector track
#
# Draws a trapezoid that fans out from the zoom region (narrow at top)
# to the full panel width (wide at bottom), visually bridging the
# Manhattan track above and the LD track below.
#
# Because the Manhattan is faceted by chromosome with
# \code{space = "free_x"}, each panel's width is proportional to its
# data range. This function computes the zoom region's position
# relative to the full faceted layout so the trapezoid top aligns
# with the highlighted rectangle on the correct chromosome panel.
#
# @param manhattanData Data from the Manhattan ggplot object.
# @param zoomChr Character chromosome ID for the zoom region.
# @param zoomStartMbp Start of zoom region in Mbp.
# @param zoomEndMbp End of zoom region in Mbp.
# @param zoomColor Color for connector lines and fill.
#
# @return A \code{ggplot} object with \code{theme_void}.
buildZoomTrack <- function(
    manhattanData,
    zoomChr,
    zoomStartMbp,
    zoomEndMbp,
    zoomColor
) {
    chrLevels <- levels(manhattanData$Chr)
    if (is.null(chrLevels)) {
        chrLevels <- unique(as.character(manhattanData$Chr))
    }
    zoomChrStr <- as.character(zoomChr)

    chrSpans <- vapply(chrLevels, function(ch) {
        chPos <- manhattanData$pos_mbp[as.character(manhattanData$Chr) == ch]
        chPos <- chPos[!is.na(chPos)]
        if (length(chPos) < 2) return(1)
        diff(range(chPos))
    }, numeric(1))

    totalSpan <- sum(chrSpans)
    if (totalSpan <= 0) totalSpan <- 1

    chrFracs  <- chrSpans / totalSpan
    cumStarts <- c(0, cumsum(chrFracs[-length(chrFracs)]))
    names(cumStarts) <- chrLevels

    chrOffset <- cumStarts[zoomChrStr]
    chrFrac   <- chrFracs[zoomChrStr]

    chrPos   <- manhattanData$pos_mbp[as.character(manhattanData$Chr) == zoomChrStr]
    chrPos   <- chrPos[!is.na(chrPos)]
    chrRange <- range(chrPos)
    chrDataSpan <- diff(chrRange)
    if (chrDataSpan <= 0) chrDataSpan <- 1

    localLeft  <- (zoomStartMbp - chrRange[1]) / chrDataSpan
    localRight <- (zoomEndMbp   - chrRange[1]) / chrDataSpan

    zoomLeft  <- max(0, min(1, chrOffset + localLeft  * chrFrac))
    zoomRight <- max(0, min(1, chrOffset + localRight * chrFrac))

    trapDF <- data.frame(
        x = c(zoomLeft, 0, 1, zoomRight),
        y = c(1, 0, 0, 1)
    )

    lineDF <- data.frame(
        x    = c(zoomLeft, zoomRight),
        xend = c(0, 1),
        y    = c(1, 1),
        yend = c(0, 0)
    )

    ggplot2::ggplot() +
        ggplot2::geom_polygon(
            data    = trapDF,
            ggplot2::aes(x = .data$x, y = .data$y),
            fill    = zoomColor,
            alpha   = 0.10
        ) +
        ggplot2::geom_segment(
            data = lineDF,
            ggplot2::aes(
                x    = .data$x,    y    = .data$y,
                xend = .data$xend, yend = .data$yend
            ),
            color     = zoomColor,
            linetype  = "dashed",
            linewidth = 0.6
        ) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::theme_void()
}


## ----
#' @title Composite Manhattan + LD track plot
#'
#' @description Stacks a pre-built Manhattan plot and a pre-built LD
#'   heatmap into a single composite figure with a zoom connector track
#'   in between. The zoom region is automatically derived from the
#'   genomic coordinates in the supplied \code{\linkS4class{LDResults}}
#'   object and highlighted on the Manhattan panel.
#'
#' @param manhattanPlot A \code{ggplot} object returned by
#'   \code{\link{plotManhattan}}.
#' @param ldPlot A \code{ggplot} object returned by
#'   \code{\link{plotLD}}.
#' @param ldObj An object of class \code{\linkS4class{LDResults}} as
#'   returned by \code{\link{linkageDiseq}}. The chromosome and
#'   position range are extracted from its \code{results} slot to
#'   determine the zoom region highlighted on the Manhattan plot.
#' @param zoomColor Color used for the highlight rectangle on the
#'   Manhattan plot and the zoom connector lines/fill. Defaults to
#'   \code{"#db0003"} (red).
#' @param heights A numeric vector of length 3 specifying the relative
#'   heights of the Manhattan, zoom connector, and LD panels.
#'   Defaults to \code{c(3, 1, 4)}.
#' @param verbose Logical. If \code{TRUE}, informational messages are
#'   printed. Defaults to \code{TRUE}.
#'
#' @return A \code{patchwork} object combining the three panels.
#'
#' @seealso \code{\link{plotManhattan}}, \code{\link{plotLD}},
#'   \code{\linkS4class{LDResults}}
#'
#' @importFrom methods is
#' @importFrom rlang .data
#'
#' @export
plotTracks <- function(
    manhattanPlot,
    ldPlot,
    ldObj,
    zoomColor = "#db0003",
    heights   = c(3, 1, 4),
    verbose   = TRUE
) {
    ## -- Validation ------------------------------------------------------
    if (!inherits(manhattanPlot, "ggplot")) {
        stop("'manhattanPlot' must be a ggplot object (e.g. from plotManhattan())")
    }
    if (!inherits(ldPlot, "ggplot")) {
        stop("'ldPlot' must be a ggplot object (e.g. from plotLD())")
    }
    if (!methods::is(ldObj, "LDResults")) {
        stop("'ldObj' must be an object of class \"LDResults\"")
    }
    if (!is.numeric(heights) || length(heights) != 3) {
        stop("'heights' must be a numeric vector of length 3")
    }

    ## -- Extract zoom region from LDResults -------------------------------
    ldDF <- as.data.frame(ldObj@results)
    loci <- unique(c(ldDF$Locus1, ldDF$Locus2))

    if (length(loci) > 1) {
        lociFreq <- table(c(ldDF$Locus1, ldDF$Locus2))
        zoomChr  <- names(which.max(lociFreq))
        if (verbose) {
            message(
                "LD data spans ", length(loci), " chromosomes; using '",
                zoomChr, "' (most frequent) for the zoom region."
            )
        }
        positions <- c(
            as.numeric(ldDF$Position1[ldDF$Locus1 == zoomChr]),
            as.numeric(ldDF$Position2[ldDF$Locus2 == zoomChr])
        )
    } else {
        zoomChr   <- loci[1]
        positions <- as.numeric(c(ldDF$Position1, ldDF$Position2))
    }

    zoomStartMbp <- min(positions) / 1e6
    zoomEndMbp   <- max(positions) / 1e6

    if (verbose) {
        message(
            "Zoom region: Chr ", zoomChr, " [",
            round(zoomStartMbp, 3), " \u2013 ",
            round(zoomEndMbp, 3), " Mbp]"
        )
    }

    ## -- Overlay highlight rectangle on Manhattan -------------------------
    chrLevels <- levels(manhattanPlot$data$Chr)
    if (is.null(chrLevels)) {
        chrLevels <- unique(as.character(manhattanPlot$data$Chr))
    }

    if (!as.character(zoomChr) %in% chrLevels) {
        warning(
            "Zoom chromosome '", zoomChr,
            "' not found in Manhattan data; ",
            "highlight rectangle will not appear."
        )
    }

    rectDF <- data.frame(
        Chr  = factor(as.character(zoomChr), levels = chrLevels),
        xmin = zoomStartMbp,
        xmax = zoomEndMbp
    )

    modifiedManhattan <- manhattanPlot +
        ggplot2::geom_rect(
            data        = rectDF,
            ggplot2::aes(
                xmin = .data$xmin, xmax = .data$xmax,
                ymin = -Inf,       ymax = Inf
            ),
            inherit.aes = FALSE,
            fill        = zoomColor,
            alpha       = 0.15,
            color       = zoomColor,
            linetype    = "dashed",
            linewidth   = 0.5
        ) +
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank()
        )

    ## -- Build zoom connector track ----------------------------------------
    zoomTrack <- buildZoomTrack(
        manhattanData = manhattanPlot$data,
        zoomChr       = zoomChr,
        zoomStartMbp  = zoomStartMbp,
        zoomEndMbp    = zoomEndMbp,
        zoomColor     = zoomColor
    )

    ## -- Assemble with patchwork ------------------------------------------
    patchwork::wrap_plots(
        modifiedManhattan,
        zoomTrack,
        ldPlot,
        ncol    = 1,
        heights = heights
    )
}
