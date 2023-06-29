## ----
#' @title Create Manhattan plot from rTASSEL association output
#'
#' @description This function allows for quick generation of a Manhattan
#'    plot from rTASSEL association statistical output data.
#'
#' @param assocRes An object of type \code{AssociationResults}
#'
#' @export
plotManhattan <- function(
        assocRes,
        trait = NULL,
        threshold = NULL,
        colors = c("#91baff", "#3e619b"),
        pltTheme = c("default", "classic"),
        showSigMarkers = FALSE
) {
    if (!is(assocRes, "AssociationResults")) {
        stop("'assocRes' parameter is not an 'AssociationResults' object")
    }

    assocType <- associationType(assocRes)
    assocStats <- switch (
        assocType,
        "GLM" = tableReport(assocRes, "GLM_Stats"),
        "MLM" = tableReport(assocRes, "MLM_Stats"),
        "FastAssoc" = tableReport(assocRes, "FastAssociation"),
        "default" = NULL
    )

    if (is.null(assocStats)) {
        stop("Association Type not defined")
    }

    plotManhattanCore(assocStats, trait, threshold, colors, pltTheme, showSigMarkers)
}


## ----
# Core visual engine
plotManhattanCore <- function(
    assocStats,
    trait,
    threshold,
    colors,
    pltTheme,
    showSigMarkers
) {
    # ## DEBUG
    # assocStats <- tableReport(tasMLM, "MLM_Stats")
    # trait <- NULL
    # threshold <- NULL
    # colors <- c("red", "black")
    # pltTheme <- "default"
    # showSigMarkers <- FALSE

    ## Coerce data frame from `DataFrame` object
    assocStats <- as.data.frame(assocStats)

    ## Filter missing chrom data
    assocStats <- assocStats[assocStats$Chr != "", ]

    ## Get unique chromosome IDs and set factor levels
    uniqChr <- unique(assocStats$Chr)
    uniqChr <- uniqChr[order(nchar(uniqChr), uniqChr)]
    assocStats$Chr <- factor(assocStats$Chr, levels = uniqChr)

    ## Filter data
    if (!is.null(trait)) {
        filtStats <- assocStats[assocStats$Trait %in% trait, ]
        filtStats$Trait <- as.factor(filtStats$Trait)
    } else {
        filtStats <- assocStats
    }


    ## Threshold highlighting
    if (!is.null(threshold)) {
        filtStats$highlight_flag <- ifelse(-log10(filtStats$p) >= threshold, TRUE, FALSE)
    } else {
        filtStats$highlight_flag <- FALSE
    }

    ## Color data (in-house use only)
    col <- c(
        "col1" = "#333333",
        "col2" = "#D34747"
    )

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

    ## Dynamic facet generation
    if (is.null(trait)) {
        trait_facet <- ggplot2::facet_grid(Trait ~ Chr, scales = "free", space = "free_x")
        main_title <- NULL
    } else if (length(trait) == 1){
        trait_facet <- ggplot2::facet_grid(. ~ Chr, scales = "free_x", space = "free_x")
        main_title <- ggplot2::ggtitle(label = paste("Trait:", paste(trait, collapse = ", ")))
    } else {
        trait_facet <- ggplot2::facet_grid(Trait ~ Chr, scales = "free", space = "free_x")
        main_title <- ggplot2::ggtitle(label = paste("Trait:", paste(trait, collapse = ", ")))
    }

    ## Transform scales
    filtStats$pos_mbp <- as.numeric(filtStats$Pos) / 1e6

    ## Plot components
    p <- ggplot2::ggplot(data = filtStats) +
        ggplot2::aes(x = .data$pos_mbp, y = -log10(.data$p)) +
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
        sig_markers +
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
# aesthetic module - "classic"
plotManhattanCoreClassic <- function() {

}


## ----
# aesthetic module - "default"
plotManhattanCoreDefault <- function() {

}


