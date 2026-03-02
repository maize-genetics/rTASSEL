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
#' @param ldBlocks An optional \code{GRanges} object specifying genomic
#'   regions to highlight as LD blocks on the plot. Each range defines a
#'   chromosome and start/end positions; sites falling within each range
#'   are outlined with a triangular border. An optional \code{label}
#'   metadata column (in \code{mcols}) adds text annotations to each
#'   block. Defaults to \code{NULL}.
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
#' @importFrom GenomicRanges seqnames start end
#' @importFrom rlang .data
#'
#' @export
plotLD <- function(tasObj,
                   ldType = c("All", "SlidingWindow"),
                   windowSize = NULL,
                   hetCalls = c("missing", "ignore", "third"),
                   plotVal = c("r2", "DPrime", "pDiseq"),
                   ldBlocks = NULL,
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

    # Sort sites by chromosome then position numerically
    sites <- unique(data.frame(
        coord = c(ldDF$coord1, ldDF$coord2),
        locus = c(ldDF$Locus1, ldDF$Locus2),
        pos   = as.numeric(c(ldDF$Position1, ldDF$Position2)),
        stringsAsFactors = FALSE
    ))
    locus_as_num <- suppressWarnings(as.numeric(sites$locus))
    if (all(!is.na(locus_as_num))) {
        sites <- sites[order(locus_as_num, sites$pos), ]
    } else {
        sites <- sites[order(sites$locus, sites$pos), ]
    }
    rownames(sites) <- NULL
    ids <- sites$coord

    # Reorder rows to match lower-triangle traversal in genomic order
    site_rank <- setNames(seq_len(nrow(sites)), sites$coord)
    r1 <- site_rank[ldDF$coord1]
    r2 <- site_rank[ldDF$coord2]
    swap <- r1 < r2
    if (any(swap)) {
        tmp <- ldDF$coord1[swap]
        ldDF$coord1[swap] <- ldDF$coord2[swap]
        ldDF$coord2[swap] <- tmp
    }
    ldDF <- ldDF[order(pmax(r1, r2), pmin(r1, r2)), ]

    ldSub    <- ldDF[, c("coord1", "coord2", plotVal)]
    ldSub    <- as.data.frame(ldSub)
    ldSubRot <- ldCellRotater(ldSub, angle)

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
        max_nchar * text_size * 0.2

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
    if (!is.null(ldBlocks)) {
        if (!inherits(ldBlocks, "GRanges")) {
            stop("'ldBlocks' must be a 'GRanges' object")
        }

        block_mcols <- S4Vectors::mcols(ldBlocks)
        has_labels  <- "label" %in% colnames(block_mcols)
        block_lw    <- max(0.5, border_lw * 2)

        for (b in seq_len(length(ldBlocks))) {
            block_chr   <- as.character(GenomicRanges::seqnames(ldBlocks)[b])
            block_start <- GenomicRanges::start(ldBlocks)[b]
            block_end   <- GenomicRanges::end(ldBlocks)[b]

            in_block <- sites$locus == block_chr &
                sites$pos >= block_start &
                sites$pos <= block_end

            if (!any(in_block)) {
                warning("Block ", b, " contains no sites, skipping")
                next
            }

            block_idx <- which(in_block)
            s <- min(block_idx)
            e <- max(block_idx)

            block_corners <- rotate(
                x = c(s, e + 1, s),
                y = c(s, e + 1, e + 1),
                angle = angle
            )

            rect_xmin <- min(block_corners$x[1], block_corners$x[2])
            rect_xmax <- max(block_corners$x[1], block_corners$x[2])
            rect_ymin <- max(block_corners$y)
            rect_ymax <- upper_margin_padding

            cell_plot <- cell_plot +
                ggplot2::annotate(
                    geom = "polygon",
                    x = block_corners$x,
                    y = block_corners$y,
                    fill = NA,
                    color = "black",
                    linewidth = block_lw
                ) +
                ggplot2::annotate(
                    geom = "rect",
                    xmin = rect_xmin,
                    xmax = rect_xmax,
                    ymin = rect_ymin,
                    ymax = rect_ymax,
                    fill = NA,
                    color = "black",
                    linewidth = block_lw
                )

            if (has_labels && !is.na(block_mcols$label[b])) {
                label_x <- (block_corners$x[1] + block_corners$x[2]) / 2
                label_y <- max(block_corners$y) - 0.15

                cell_plot <- cell_plot +
                    ggplot2::annotate(
                        geom = "text",
                        x = label_x,
                        y = label_y,
                        label = block_mcols$label[b],
                        size = text_size + 0.5,
                        fontface = "bold"
                    )
            }
        }
    }

    return(cell_plot)
}
