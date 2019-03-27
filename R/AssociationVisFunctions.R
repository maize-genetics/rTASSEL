#--------------------------------------------------------------------
# Script Name:   AssociationVisFunctions.R
# Description:   Visualization functions for association analyses
# Author:        Brandon Monier
# Created:       2019-03-13 at 13:27:30
# Last Modified: 2019-03-13 at 16:20:06
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions for
#    visualizing association output of TASSEL.
#--------------------------------------------------------------------

manhattanPlot <- function(assocStats,
                          trait,
                          threshold,
                          colors = c("#91baff", "#3e619b"),
                          showSigMarkers = FALSE) {

    ## Get chromosome count
    chrs <- length(levels(assocStats$Chr))

    ## Filter data
    filtglmStats <- assocStats %>%
        dplyr::filter(Trait == trait) %>%
        dplyr::mutate(highlight_flag = ifelse(-log10(p) >= threshold, T, F))

    ## Color data
    col <- c(
        "col1" = "#333333",
        "col2"   = "#D34747"
    )

    ## Plot components
    p <- ggplot2::ggplot(data = filtglmStats) +
        ggplot2::aes(x = Pos, y = -log10(p)) +
        ggplot2::geom_point(size = 0.8, na.rm = TRUE) +
        ggplot2::aes(color = Chr) +
        ggplot2::geom_point(
            data = dplyr::filter(
                filtglmStats, highlight_flag == T
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
                        yes = as.character(Marker),
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
