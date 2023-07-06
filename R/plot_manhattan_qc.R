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
#' @param threshold User-defined \eqn{-log_{10}(p)}-value threshold for
#'    significant marker determination. Once specified any marker points
#'    higher than this line will be highlighted.
#'
#' @return Returns a \code{ggplot2} object
#'
#' @export
plotManhattanQC <- function(
    assocRes,
    gr,
    trait = NULL,
    window = NULL,
    threshold = NULL
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
        "assocStats" = assocStats,
        "gr"         = gr,
        "threshold"  = threshold,
        "trait"      = trait,
        "window"     = window
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
    filtStats <- filterByGRanges(filtStats, gr, window)

    ## Dynamic facet generation
    if (is.null(trait) || length(trait) > 1) {
        trait_facet <- ggplot2::facet_grid(Trait ~ range_id, scales = "free")
        main_title <- NULL
    } else {
        trait_facet <- ggplot2::facet_grid(. ~ range_id, scales = "free_x")
        main_title <- ggplot2::ggtitle(label = paste("Trait:", paste(trait, collapse = ", ")))
    }

    ## Plot components
    p <- ggplot2::ggplot(data = filtStats) +
        ggplot2::aes(x = .data$pos_mbp, y = -log10(.data$p)) +
        ggplot2::geom_point(size = 0.8, na.rm = TRUE) +
        # ggplot2::aes(color = .data$Chr) +
        # ggplot2::scale_color_manual(
        #     values = rep(colors, length(levels(filtStats$Chr)))
        # ) +
        # ggplot2::geom_vline(xintercept = 5000000) +
        # ggplot2::geom_vline(xintercept = 6000000) +
        ggplot2::geom_hline(yintercept = threshold, linetype = "dashed") +
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
filterByGRanges <- function(assocStats, gr, window) {
    ## Coerce to data frame (if GRanges object)
    grDf <- as.data.frame(gr)

    ## Sanity check for columns
    neededCols <- c("seqnames", "start", "end", "range_id")
    for (col in neededCols) {
        assocStatsColumnChecker(col, grDf, neededCols)
    }

    ## Drop possible factors
    grDf$seqnames <- as.character(grDf$seqnames)

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
            if ("canonical_id" %in% colnames(tmpRow)) {
                filtStats$range_id <- tmpRow$canonical_id
            } else {
                filtStats$range_id <- tmpRow$range_id
            }

            filtStats$range_start <- tmpRow$start
            filtStats$range_end   <- tmpRow$end
            dfLs[[row]] <- filtStats
        } else {
            message("No association results found for: ", tmpRow$range_id)
        }
    }

    dfLs <- do.call("rbind", dfLs)
    if (nrow(dfLs) > 0) {
        return(dfLs)
    } else {
        stop("No association results found for any trait")
    }
}















