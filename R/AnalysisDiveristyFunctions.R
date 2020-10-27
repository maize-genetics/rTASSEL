#--------------------------------------------------------------------
# Script Name:   AnalysisDiversityFunctions.R
# Description:   Functions for caclulating linkage disequilibrium
# Author:        Brandon Monier
# Created:       2019-06-27 at 10:00:00
# Last Modified: 2020-06-23 at 11:50:36
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#   The main purpose of this Rscript to host functions necessary for
#   methods pertaining to calculating LD from rTASSEL genotype
#   datasets.
#--------------------------------------------------------------------

#' @title Calculate linkage disequilibrium from an rTASSEL genotype
#'   dataset.
#'
#' @description Generates a linkage disequilibrium (LD) data set from SNP data.
#'
#' @name linkageDiseq
#' @rdname linkageDiseq
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype} that
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
#' @return Returns a \code{DataFrame}-based data frame.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava J
#' @importFrom S4Vectors DataFrame
#'
#' @export
linkageDiseq <- function(tasObj,
                         ldType = c("SlidingWindow", "All"),
                         windowSize = NULL,
                         hetCalls = c("missing", "ignore", "third"),
                         verbose = TRUE) {

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

    tableReportToDF(ldObj)
}


