## ----
#' @title Create a QC Manhattan plots from rTASSEL association output
#'
#' @description This function allows for quick generation of a QC
#'    plot from rTASSEL association statistical output data. The main goal
#'    of this function is to provide "zoomed" in Manhattan plots that
#'    typically fall within genomic ranges of interest plus a flanking
#'    window up- and downstream of the ranges.
#'
#' @param assocRes An object of type \code{AssociationResults}
#' @param gr Genomic ranges of interest. Can be passed as a \code{GRanges}
#'    object or a \code{data.frame} object.
#' @param trait Which phenotypic trait do you want to plot? If set to
#'    \code{NULL}, this will generate a faceted plot with all mapped traits.
#' @param window Window size (base-pairs) flanking surround reference range.
#'    Defaults to \code{100000} (100,000 base-pairs).
#' @param threshold User-defined \eqn{-log_{10}(p)}-value threshold for
#'    significant marker determination. Once specified any marker points
#'    higher than this line will be highlighted.
#' @param classicNames Do you want to plot classical gene names instead?
#'    NOTE: this will need a \code{classical_id} column for the "ranges of
#'    interest" data. Defaults to \code{FALSE}.
#' @param verbose Should messages be printed to console? Defaults to
#'    \code{FALSE}.
#'
#' @return Returns a \code{ggplot2} object
#'
#' @export
plotManhattanQC <- function(
    assocRes,
    trait        = NULL,
    gr,
    window       = 100000,
    threshold    = NULL,
    classicNames = FALSE,
    verbose      = FALSE
) {
    if (!is(assocRes, "AssociationResults")) {
        stop(
            "The object '", deparse(substitute(assocRes)),
            "' is not an 'AssociationResults' object"
        )
    }

    traitValidityChecker(trait, assocRes)

    assocStats <- switch (
        associationType(assocRes),
        "GLM"       = tableReport(assocRes, "GLM_Stats"),
        "MLM"       = tableReport(assocRes, "MLM_Stats"),
        "FastAssoc" = tableReport(assocRes, "FastAssociation"),
        "default"   = NULL
    )

    if (is.null(assocStats)) {
        stop("Association Type not defined")
    }

    # Make a more robust/'OOP' object to pass
    pManQCCoreParams <- list(
        "assocStats"   = assocStats,
        "classicNames" = classicNames,
        "gr"           = gr,
        "threshold"    = threshold,
        "trait"        = trait,
        "window"       = window,
        "verbose"      = verbose
    )

    plotManhattanQCCore(pManQCCoreParams)
}


## ----
#' @title Core visual engine for QQ plotting
#' @param params A list of parameter variables
#' @importFrom rlang .data
plotManhattanQCCore <- function(params) {

    ## Parse parameters
    assocStats <- params$assocStats
    gr         <- params$gr
    threshold  <- params$threshold
    trait      <- params$trait
    window     <- params$window

    ## "Prime" and filter data for plotting
    filtStats <- primeManhattanData(params)

    ## Filter based on ranges of interest
    filtStats <- filterByGRanges(filtStats, params)

    ## Convert ranges to mapping data
    grDfAes <- grToAesMappingDf(params)

    ## Dynamic facet generation
    if (is.null(trait) || length(trait) > 1) {
        trait_facet <- ggplot2::facet_grid(Trait ~ range_id, scales = "free")
        main_title  <- NULL
    } else {
        trait_facet <- ggplot2::facet_grid(. ~ range_id, scales = "free_x")
        main_title  <- ggplot2::ggtitle(label = paste("Trait:", paste(trait, collapse = ", ")))
    }

    ## Plot components
    p <- ggplot2::ggplot() +
        ggplot2::geom_point(
            data = filtStats,
            mapping = ggplot2::aes(
                x = .data$pos_mbp,
                y = -log10(.data$p)
            ),
            size = 0.8,
            na.rm = TRUE
        ) +
        ggplot2::geom_vline(
            data = grDfAes$aes_vlin,
            mapping = ggplot2::aes(
                xintercept = .data$x_pos
            ),
            linetype = "dashed"
        ) +
        ggplot2::geom_hline(yintercept = threshold, linetype = "dashed") +
        ggplot2::geom_rect(
            data    = grDfAes$aes_rect,
            fill    = "blue",
            alpha   = 0.2,
            mapping = ggplot2::aes(
                xmin  = .data$xmin,
                xmax  = .data$xmax,
                ymin  = .data$ymin,
                ymax  = .data$ymax
            )
        ) +
        # Hack 1 - "set" lower x bound
        ggplot2::geom_blank(
            data = grDfAes$aes_vlin,
            mapping = ggplot2::aes(
                y = 0,
                x = .data$x_start_bound
            )
        ) +
        # Hack 2 - "set" upper x bound
        ggplot2::geom_blank(
            data = grDfAes$aes_vlin,
            mapping = ggplot2::aes(
                y = 0,
                x = .data$x_end_bound
            )
        ) +
        ggplot2::xlab("SNP Position (Mbp)") +
        ggplot2::ylab(bquote(~-log[10]~ '('*italic(p)*'-value)')) +
        main_title +
        ggplot2::theme_bw() +
        trait_facet +
        ggplot2::theme(
            legend.position    = "none",
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank()
        )

    return(p)
}


## ----
filterByGRanges <- function(assocStats, params) {
    ## Parse parameters
    gr           <- params$gr
    window       <- params$window
    classicNames <- params$classicNames
    verbose      <- params$verbose

    ## Coerce to data frame (if GRanges object)
    grDf <- as.data.frame(gr)

    ## Sanity check for columns
    neededCols <- c("seqnames", "start", "end", "range_id")
    checkForValidColumns(grDf, neededCols)

    ## Drop possible factors
    grDf$seqnames <- as.character(grDf$seqnames)

    ## Check for classic ID
    idType <- checkForClassicId(classicNames, grDf, muteWarnings = FALSE)

    ## Convoluted base R filtering!
    dfLs <- vector("list", length = nrow(grDf))
    for (row in seq_len(nrow(grDf))) {
        tmpRow <- grDf[row, ]
        if ((tmpRow$start - window) < 1) {
            windowStart <- 1
            windowEnd   <- tmpRow$end + window
        } else {
            windowStart <- tmpRow$start - window
            windowEnd   <- tmpRow$end + window
        }
        filtStats <- assocStats[
            assocStats$Chr == tmpRow$seqnames &
            assocStats$Pos >= windowStart &
            assocStats$Pos <= windowEnd,
        ]
        if (nrow(filtStats) != 0) {
            filtStats$range_id    <- tmpRow[, idType]
            filtStats$range_start <- tmpRow$start
            filtStats$range_end   <- tmpRow$end
            dfLs[[row]] <- filtStats
        } else {
            if (verbose) message("No association results found for: ", tmpRow[, idType])
        }
    }

    dfLs <- do.call("rbind", dfLs)
    if (!is.null(dfLs)) {
        return(dfLs)
    } else {
        stop("No association results found for any trait")
    }
}


## ----
checkForClassicId <- function(
    classicNames,
    grDf,
    muteWarnings = TRUE
) {
    ID_TYPE <- list(
        "DEFAULT" = "range_id",
        "CLASSIC" = "classical_id"
    )

    if (classicNames) {
        if (!"classical_id" %in% colnames(grDf)) {
            if (!muteWarnings) {
                warning("'classical_id' parameter not found. Defaulting to 'range_id'")
            }
            idType <- ID_TYPE$DEFAULT
        } else {
            idType <- ID_TYPE$CLASSIC
        }
    } else {
        idType <- ID_TYPE$DEFAULT
    }

    return(idType)
}


## ----
grToAesMappingDf <- function(params) {
    ## Parse parameters
    gr           <- params$gr
    window       <- params$window
    classicNames <- params$classicNames
    scale        <- 1e6

    grDf <- as.data.frame(gr)

    ## Sanity check for columns
    neededCols <- c("seqnames", "start", "end", "range_id")
    checkForValidColumns(grDf, neededCols)

    idType <- checkForClassicId(classicNames, grDf)

    aesDf <- data.frame(
        xmin     = grDf$start / scale,
        xmax     = grDf$end / scale,
        ymin     = -Inf,
        ymax     = Inf,
        range_id = grDf[, idType]
    )

    aesVLineDf <- data.frame(
        x_pos         = ((grDf$start / scale) + (grDf$end / scale)) * 0.5,
        x_start_bound = (grDf$start / scale) - (window / scale),
        x_end_bound   = (grDf$end / scale) + (window / scale),
        range_id      = grDf[, idType]
    )

    return(
        list(
            "aes_rect" = aesDf,
            "aes_vlin" = aesVLineDf
        )
    )
}


