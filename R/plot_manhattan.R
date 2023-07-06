## ----
#' @title Create Manhattan plot from rTASSEL association output
#'
#' @description This function allows for quick generation of a Manhattan
#'    plot from rTASSEL association statistical output data.
#'
#' @param assocRes An object of type \code{AssociationResults}
#' @param trait Which phenotypic trait do you want to plot? If set to
#'    \code{NULL}, this will generate a faceted plot with all mapped traits
#' @param threshold User-defined \eqn{-log_{10}(p)}-value threshold for
#'    significant marker determination. Once specified any marker points
#'    higher than this line will be highlighted.
#' @param colors A vector of \code{character} colors used for differentiating
#'    multiple chromosomes. Defaults to 2 shades of blue.
#' @param pltTheme A given plot
#'
#' @return Returns a \code{ggplot2} object
#'
#' @export
plotManhattan <- function(
        assocRes,
        trait = NULL,
        threshold = NULL,
        colors = c("#91baff", "#3e619b"),
        pltTheme = c("default", "classic")
) {
    if (!is(assocRes, "AssociationResults")) {
        stop(
            "The object '", deparse(substitute(assocRes)),
            "' is not an 'AssociationResults' object"
        )
    }

    if (!is.null(trait) && !trait %in% traitNames(assocRes)) {
        stop("Trait: '", trait, "' is not a valid trait ID found in this object")
    }

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

    plotManhattanCore(assocStats, trait, threshold, colors, pltTheme)
}


## ----
#' @title Core visual engine
#' @importFrom rlang .data
plotManhattanCore <- function(
    assocStats,
    trait,
    threshold,
    colors,
    pltTheme
) {

    filtStats <- primeManhattanData(assocStats, trait, threshold)

    ## Dynamic facet generation
    if (is.null(trait) || length(trait) > 1) {
        trait_facet <- ggplot2::facet_grid(Trait ~ Chr, scales = "free", space = "free_x")
        main_title <- NULL
    } else {
        trait_facet <- ggplot2::facet_grid(. ~ Chr, scales = "free_x", space = "free_x")
        main_title <- ggplot2::ggtitle(label = paste("Trait:", paste(trait, collapse = ", ")))
    }

    ## Plot components
    p <- ggplot2::ggplot(data = filtStats) +
        ggplot2::aes(x = .data$pos_mbp, y = -log10(.data$p)) +
        ggplot2::geom_point(size = 0.8, na.rm = TRUE) +
        ggplot2::aes(color = .data$Chr) +
        ggplot2::scale_color_manual(
            values = rep(colors, length(levels(filtStats$Chr)))
        ) +
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
# plotManhattan dataframe primer
primeManhattanData <- function(
    assocStats,
    trait,
    threshold
) {
    ## Sanity check coerce data frame from assocStats object
    # NOTE - probably don't need to do this right now but a security
    #        feature if I change TableReport classes in the future...
    assocStats <- as.data.frame(assocStats)

    ## Sanity check for columns
    neededCols <- c("Chr", "Pos", "Trait")
    for (col in neededCols) {
        assocStatsColumnChecker(col, assocStats, neededCols)
    }

    ## Coerce Chr and Pos vectors
    assocStats$Chr <- as.character(assocStats$Chr)
    assocStats$Pos <- as.numeric(assocStats$Pos)

    ## Filter missing data
    assocStats <- assocStats[assocStats$Chr != "", ]
    assocStats <- assocStats[assocStats$Trait != "", ]

    ## Get unique chromosome IDs and set factor levels
    uniqChr <- unique(assocStats$Chr)
    uniqChr <- uniqChr[order(nchar(uniqChr), uniqChr)]
    assocStats$Chr <- factor(assocStats$Chr, levels = uniqChr)

    ## Transform scales
    assocStats$pos_mbp <- as.numeric(assocStats$Pos) / 1e6

    ## Filter data
    if (!is.null(trait)) {
        filtStats <- assocStats[assocStats$Trait %in% trait, ]
        filtStats$Trait <- as.factor(filtStats$Trait)
    } else {
        filtStats <- assocStats
    }

    ## Threshold highlighting
    if (!is.null(threshold)) {
        filtStats$highlight_flag <- ifelse(
            test = -log10(filtStats$p) >= threshold,
            yes  = TRUE,
            no   = FALSE
        )
    } else {
        filtStats$highlight_flag <- FALSE
    }

    return(filtStats)
}


## ----
# sanity checker for association results
assocStatsColumnChecker <- function(x, assocStats, neededCols) {
    if (!x %in% colnames(assocStats)) {
        stop(
            "'", x, "' column not found in stats dataframe - need at least: ",
            paste(neededCols, collapse = ",")
        )
    }
}


