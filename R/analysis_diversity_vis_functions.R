## ----
#' @title Linkage desequilibrium visualization application
#'
#' @description Calculates linkage disequilibrium (LD) and runs an
#'   interactive Java visualizer for LD results.
#'
#' @name ldJavaApp
#' @rdname ldJavaApp
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
#'   current site.
#' @param windowSize What size do you want your LD analysis window? If you
#'   have chosen \code{SlidingWindow} for the \code{ldType} parameter, you
#'   will need to specify window size.
#' @param hetCalls How should heterozygous calls be handled? Current options
#'   are \code{"ignore"} (ignore altogether), \code{"missing"}
#'   (set to missing), and \code{"third"} (treat as third state).
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
#' @seealso \code{\link{linkageDiseq}}, \code{\link{ldPlot}}
#'
#' @return Returns a Java-based visualization application.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava J
#' @importFrom rJava .jnew
#'
#' @export
# nocov start
ldJavaApp <- function(tasObj,
                      ldType = c("SlidingWindow", "All"),
                      windowSize = NULL,
                      hetCalls = c("missing", "ignore", "third"),
                      verbose = TRUE) {
    warnMsg <- paste0("The function 'ldJavaApp()' will be deprecated soon.")
    message(warnMsg)

    # Logic - Check for TasselGenotypePhenotype class
    if (!is(tasObj, "TasselGenotypePhenotype")) {
        stop("tasObj is not of class \"TasselGenotypePhenotype\"")
    }

    # Logic - Check to see if TASSEL object has a genotype table
    if (rJava::is.jnull(tasObj@jGenotypeTable)) {
        stop("tasObj does contain a Genotype object")
    }

    # Logic - Check for available parameters
    hetCalls <- match.arg(hetCalls)
    ldType <- match.arg(ldType)

    # Logic - Add warning for all
    if (ldType == "All") {
        if (verbose) message("This *might* produce a massive matrix. You have been warned.")
        windowSize <- -1
    }

    # Logic - check for ldType and windowSize compatability
    if (ldType == "SlidingWindow" & is.null(windowSize)) {
        if (verbose) message("`windowSize` is not set - setting to `1`")
        windowSize <- 1
    }

    # Get TASSEL generate R code plugin
    jRC <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")

    # Run LD
    if (verbose) message("Calculating LD...")
    ldObj <- jRC$linkageDiseq(
        tasObj@jGenotypeTable,  # TASSEL genotype table
        ldType,                 # LD type parameter
        as.integer(windowSize), # Window size
        hetCalls                # heterozygous calls
    )

    ## Call LD dialog ----
    ldPlug <- rJava::.jnew(
        class = "net/maizegenetics/analysis/popgen/LinkageDiseqDisplayPlugin",
        rJava::.jnew("java/awt/Frame"),
        TRUE
    )
    ldDialog <- rJava::.jnew(
        class = "net/maizegenetics/analysis/popgen/LinkageDiseqDisplayDialog",
        ldPlug,
        ldObj
    )
    ldDialog$show()
}
# nocov end



#' @title Linkage disequilibrium plot
#'
#' @description Calculates linkage disequilibrium (LD) and generates
#'   a static plot using \code{ggplot2} graphics.
#'
#' @name ldPlot
#' @rdname ldPlot
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
#' @seealso \code{\link{linkageDiseq}}, \code{\link{ldJavaApp}}
#'
#' @return Returns a \code{ggplot2} object.
#'
#' @importFrom rlang .data
#'
#' @export
ldPlot <- function(tasObj,
                   ldType = c("All", "SlidingWindow"),
                   windowSize = NULL,
                   hetCalls = c("missing", "ignore", "third"),
                   plotVal = c("r2", "DPrime", "pDiseq"),
                   verbose = TRUE) {
    # House-keeping parameters
    angle                <- 135
    text_scale_buffer    <- 1.5
    upper_margin_padding <- 3

    # Plot logic handling
    plotVal <- match.arg(plotVal)
    ldType  <- match.arg(ldType)

    # Generate initial data frame
    ldDF <- linkageDiseq(
        tasObj     = tasObj,
        ldType     = ldType,
        windowSize = windowSize,
        hetCalls   = hetCalls,
        verbose    = verbose
    )

    # Subset LD data frame
    if (plotVal == "r2") plotVal <- "R^2"
    ldDF$coord1 <- paste0(ldDF$Locus1, "_", ldDF$Position1)
    ldDF$coord2 <- paste0(ldDF$Locus2, "_", ldDF$Position2)
    ldSub        <- ldDF[, c("coord1", "coord2", plotVal)]
    ldSub        <- as.data.frame(ldSub)

    # Rotate coordinates for plotting
    ldSubRot <- ldCellRotater(ldSub, angle)

    # Get IDs and rotate to new position
    ids <- unique(c(ldSub$coord1, ldSub$coord2))
    ids <- ids[order(ids)]
    id_coord <- list(
        "x" = seq(1.5, 1.5 + length(ids) - 1, 1),
        "y" = seq(0.5, 0.5 + length(ids) - 1, 1)
    )
    id_coord <- rotate(id_coord$x, id_coord$y, angle)

    # Generate "clean" legend label
    legend_lab <- switch(
        EXPR = plotVal,
        "R^2" = bquote(italic(r)^2),
        "DPrime" = bquote(italic(D)*"'"),
        "pDiseq" = bquote(italic(p)*'-value')
    )

    # Make ggplot2 object
    cell_plot <- ggplot2::ggplot(data = ldSubRot) +
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$val, group = .data$group) +
        ggplot2::geom_polygon(color = "white", linewidth = 1) +
        ggplot2::annotate(
            geom = "text",
            x = id_coord$x,
            y = id_coord$y + (abs(id_coord$y) * text_scale_buffer),
            label = ids,
            angle = 90
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


