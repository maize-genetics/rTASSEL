#--------------------------------------------------------------------
# Script Name:   AssociationVisFunctions.R
# Description:   Visualization functions for association analyses
# Author:        Brandon Monier
# Created:       2019-03-13 at 13:27:30
# Last Modified: 2020-06-15 at 15:27:18
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions for
#    visualizing association output of TASSEL.
#--------------------------------------------------------------------

#' @title Create Manhattan plot from rTASSEL association output
#'
#' @description This function allows for quick generation of a Manhattan
#'    plot from rTASSEL association statistical output data.
#'
#' @name manhattanPlot
#' @rdname manhattanPlot
#'
#' @param assocStats rTASSEL association statistical output. This output is
#'    derived from association output dat frames ending in \code{_Stats}.
#'    This data contains the necessary column variables for this plot to work.
#' @param trait Which phenotypic trait do you want to plot?
#' @param threshold User-defined \eqn{-log_{10}(p)}-value threshold for significant
#'    marker determination. Once specified any marker points higher than this
#'    line will be highlighted.
#' @param colors A vector of colors used for differentiating multiple
#'    chromosomes. Defaults to 2 shades of blue.
#' @param showSigMarkers Should significantly highlighted markers be labelled?
#'    If \code{TRUE}, marker points will labelled with their respective
#'    marker labels. Defaults to \code{FALSE}.
#' @param showRug Should a rug plot be plotted under the Manhattan
#'    plot? Defaults to \code{FALSE}.
#'
#' @return Returns a \code{ggplot2} object
#'
#' @importFrom rlang .data
#'
#' @export
manhattanPlot <- function(assocStats,
                          trait,
                          threshold,
                          colors = c("#91baff", "#3e619b"),
                          showSigMarkers = FALSE,
                          showRug = FALSE) {

    ## Coerce data frame from `DataFrame` object
    assocStats <- as.data.frame(assocStats)

    ## Logic check for threshold
    if (missing(threshold)) stop("Please enter a numeric threshold value.")

    ## Get unique chromosome IDs and set factor levels
    uniqChr <- unique(assocStats$Chr)
    uniqChr <- uniqChr[order(nchar(uniqChr), uniqChr)]
    assocStats$Chr <- factor(assocStats$Chr, levels = uniqChr)

    ## Filter data
    filtStats <- assocStats[which(assocStats$Trait == trait), ]
    filtStats$highlight_flag <- ifelse(-log10(filtStats$p) >= threshold, TRUE, FALSE)

    ## Color data (in-house use only)
    col <- c(
        "col1" = "#333333",
        "col2" = "#D34747"
    )

    ## ggplot2 component logic - rug plots
    if (!showRug) {
        rug_comp <- NULL
    } else {
        rug_comp <- ggplot2::geom_rug(
            sides  = "b",
            colour = "black",
            alpha  = 0.1
        )
    }

    ## ggplot2 component logic - show sig. markers
    if (!showSigMarkers) {
        sig_markers <- NULL
    } else {
        sig_markers <- ggplot2::geom_text(
            ggplot2::aes(
                label = ifelse(
                    test = -log10(.data$p) >= threshold,
                    yes  = as.character(.data$Marker),
                    no   = ""
                )
            ),
            color = col[["col2"]],
            size  = 3,
            hjust = "inward",
            vjust = "inward",
            na.rm = TRUE
        )
    }

    ## Plot components
    p <- ggplot2::ggplot(data = filtStats) +
        ggplot2::aes(x = .data$Pos, y = -log10(.data$p)) +
        ggplot2::geom_point(size = 0.8, na.rm = TRUE) +
        ggplot2::aes(color = .data$Chr) +
        ggplot2::geom_point(
            data  = filtStats[which(filtStats$highlight_flag == TRUE), ],
            color = col[["col1"]],
            size  = 2
        ) +
        ggplot2::scale_color_manual(
            values = rep(colors, length(uniqChr))
        ) +
        rug_comp +
        sig_markers +
        ggplot2::geom_hline(yintercept = threshold, linetype = "dashed") +
        ggplot2::xlab("Position") +
        ggplot2::ylab(bquote(~-log[10]~ '('*italic(p)*'-value)')) +
        ggplot2::ggtitle(label = paste("Trait:", trait)) +
        ggplot2::facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x        = ggplot2::element_text(angle = 90, hjust = 1),
            legend.position    = "none",
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank()
        )

    return(p)
}
