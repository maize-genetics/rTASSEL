#--------------------------------------------------------------------
# Script Name:   AssociationVisFunctions.R
# Description:   Visualization functions for association analyses
# Author:        Brandon Monier
# Created:       2019-03-13 at 13:27:30
# Last Modified: 2019-04-04 at 23:32:09
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
#' @param threshold User-defined \eqn{p}-value threshold for significant
#'    marker determination. Once specified any marker points higher than this
#'    line will be highlighted.
#' @param colors A vector of colors used for differentiating multiple
#'    chromosomes. Defaults to 2 shades of blue.
#' @param showSigMarkers Should significantly highlighted markers be labelled?
#'    If \code{TRUE}, marker points will labelled with their respective
#'    marker labels. Defaults to \code{FALSE}.
#'
#' @return Returns a \code{ggplot2} object
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
manhattanPlot <- function(assocStats,
                          trait,
                          threshold,
                          colors = c("#91baff", "#3e619b"),
                          showSigMarkers = FALSE) {

    ## Get chromosome count
    chrs <- length(levels(assocStats$Chr))

    ## Filter data
    filtglmStats <- assocStats %>%
        dplyr::filter(.data$Trait == trait) %>%
        dplyr::mutate(highlight_flag = ifelse(-log10(p) >= threshold, T, F))

    ## Color data
    col <- c(
        "col1" = "#333333",
        "col2"   = "#D34747"
    )

    ## Plot components
    p <- ggplot2::ggplot(data = filtglmStats) +
        ggplot2::aes(x = .data$Pos, y = -log10(.data$p)) +
        ggplot2::geom_point(size = 0.8, na.rm = TRUE) +
        ggplot2::aes(color = .data$Chr) +
        ggplot2::geom_point(
            data = dplyr::filter(
                filtglmStats, .data$highlight_flag == T
            ),
            color = col[["col1"]],
            size = 2
        ) +
        ggplot2::scale_color_manual(
            values = rep(colors, chrs)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 90, hjust = 1
            ),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank()
        ) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::geom_hline(yintercept = threshold, linetype = "dashed") +
        ggplot2::xlab("Position") +
        ggplot2::ylab(bquote(~-log[10]~ '('*italic(p)*'-value)')) +
        ggplot2::geom_rug(sides = "b", colour = col[["black"]], alpha = 0.1) +
        ggplot2::ggtitle(label = paste("Trait:", trait)) +
        ggplot2::facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
        if (!showSigMarkers) {
            NULL
        } else {
            ggplot2::geom_text(
                ggplot2::aes(
                    label = ifelse(
                        test = -log10(p) >= threshold,
                        yes = as.character(.data$Marker),
                        no = ""
                    )
                ),
                color = col[["col2"]],
                size = 3,
                hjust = "inward",
                vjust = "inward",
                na.rm = TRUE
            )
        }
    return(p)
}
