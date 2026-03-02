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
#' @importFrom rlang .data
#'
#' @export
plotLD <- function(tasObj,
                   ldType = c("All", "SlidingWindow"),
                   windowSize = NULL,
                   hetCalls = c("missing", "ignore", "third"),
                   plotVal = c("r2", "DPrime", "pDiseq"),
                   verbose = TRUE) {
    angle     <- 135
    label_gap <- 0.2

    plotVal <- match.arg(plotVal)
    ldType  <- match.arg(ldType)

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
    ldSub        <- ldDF[, c("coord1", "coord2", plotVal)]
    ldSub        <- as.data.frame(ldSub)

    ldSubRot <- ldCellRotater(ldSub, angle)

    ids <- sort(unique(c(ldSub$coord1, ldSub$coord2)))
    n_sites <- length(ids)
    id_coord <- list(
        "x" = seq(1.5, 1.5 + n_sites - 1, 1),
        "y" = seq(0.5, 0.5 + n_sites - 1, 1)
    )
    id_coord <- rotate(id_coord$x, id_coord$y, angle)

    border_lw <- min(1, 10 / n_sites)
    text_size <- max(1.5, min(3.5, 40 / n_sites))

    max_nchar <- max(nchar(ids))
    upper_margin_padding <- max(id_coord$y) + label_gap +
        max_nchar * text_size * 0.1

    legend_lab <- switch(
        EXPR = plotVal,
        "R^2" = bquote(italic(r)^2),
        "DPrime" = bquote(italic(D)*"'"),
        "pDiseq" = bquote(italic(p)*'-value')
    )

    cell_plot <- ggplot2::ggplot(data = ldSubRot) +
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$val, group = .data$group) +
        ggplot2::geom_polygon(color = "white", linewidth = border_lw) +
        ggplot2::annotate(
            geom = "text",
            x = id_coord$x,
            y = id_coord$y + label_gap,
            label = ids,
            angle = 90,
            hjust = 0,
            size = text_size
        ) +
        ggplot2::ylim(min(ldSubRot$y), upper_margin_padding) +
        ggplot2::scale_x_reverse() +
        ggplot2::scale_fill_continuous(
            name = legend_lab,
            type = "viridis"
        ) +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(
            legend.position = "bottom"
        )
    return(cell_plot)
}
