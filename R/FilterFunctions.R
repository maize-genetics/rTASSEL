#---------------------------------------------------------------------
# Script Name:   FilterFunctions.R
# Description:   Functions for TASSEL filtration
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2020-10-26 at 13:14:55
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#   The main purpose of this Rscript to host functions necessary for
#   filtration methods of TASSEL genotype tables
#--------------------------------------------------------------------

#' @title Filter genotype table by sites
#'
#' @description This function will filter R objects of
#'    \code{TasselGenotypePhenotype} class containing genotype tables.
#'    The parameters for this function are derived from TASSEL's
#'    \code{FilterSiteBuilder} plugin.
#'
#' @name filterGenotypeTableSites
#' @rdname filterGenotypeTableSites
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#' @param siteMinCount Site minimum count of alleles not unknown. Can range
#'    from 0 to inf. Defaults to 0.
#' @param siteMinAlleleFreq Site minimum minor allele frequency. Can range
#'    from 0 to 1.0. Defaults to 0.0.
#' @param siteMaxAlleleFreq Site maximum minor allele frequency. Can range
#'    from 0 to 1.0. Defaults to 1.0.
#' @param minHeterozygous Min heterozygous proportion. Can range from 0 to 1.0.
#'    Defaults to 0.0.
#' @param maxHeterozygous Max heterozygous proportion. Can range from 0 to 1.0.
#'    Defaults to 1.0.
#' @param removeMinorSNPStates Remove minor SNP states. Defaults to
#'    \code{FALSE}.
#' @param removeSitesWithIndels Remove sites containing an indel
#'    (\code{+} or \code{-}). Defaults to \code{FALSE}.
#' @param siteRangeFilterType True if filtering by site numbers. False if
#'    filtering by chromosome and position. Options are
#'    \code{none}, \code{sites}, or \code{position}. Defaults to \code{none}.
#' @param startSite The start site. Defaults to 0.
#' @param endSite The end site. Defaults to 0.
#' @param startChr Start chromosome for site filtration range if \code{position}
#'    is chosen from \code{siteRangeFilterType}. Needs end chromosome
#'    (\code{endChr}) to work.
#' @param endChr End chromosome for site filtration range if \code{position}
#'    is chosen from \code{siteRangeFilterType}. Needs start chromosome
#'    (\code{endChr}) to work.
#' @param startPos Physical start position (bp) for filtration range if
#'    \code{position} is chosen from \code{siteRangeFilterType}. If
#'    \code{NULL}, the first physical position in the data set will be
#'    chosen.
#' @param endPos Physical end position (bp) for filtration range if
#'    \code{position} is chosen from \code{siteRangeFilterType}. If
#'    \code{NULL}, the last physical position in the data set will be
#'    chosen.
#' @param chrPosFile An optional chromosome position file path of
#'    \code{character} class. Defaults to \code{NULL}. \strong{Note:}
#'    a chromosome position file must contain correct formatting
#'    (e.g. a two column file with the header of
#'    \code{c("Chromosome", "Position")}).
#' @param bedFile An optional BED coordinate file path of
#'    \code{character} class. Defaults to \code{NULL}.
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#' @export
filterGenotypeTableSites <- function(tasObj,
                                     siteMinCount = 0,
                                     siteMinAlleleFreq = 0.0,
                                     siteMaxAlleleFreq = 1.0,
                                     minHeterozygous = 0.0,
                                     maxHeterozygous = 1.0,
                                     removeMinorSNPStates = FALSE,
                                     removeSitesWithIndels = FALSE,
                                     siteRangeFilterType = c("none", "sites", "position"),
                                     startSite = NULL,
                                     endSite = NULL,
                                     startChr = NULL,
                                     startPos = NULL,
                                     endChr = NULL,
                                     endPos = NULL,
                                     chrPosFile = NULL,
                                     bedFile = NULL) {

    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    # Range check
    if (siteMinAlleleFreq > 1 || siteMinAlleleFreq < 0) {
        stop("siteMinAlleleFreq parameter is out of range")
    }
    if (siteMaxAlleleFreq > 1 || siteMaxAlleleFreq < 0) {
        stop("siteMaxAlleleFreq parameter is out of range")
    }
    if (minHeterozygous > 1 || minHeterozygous < 0) {
        stop("minHeterozygous parameter is out of range")
    }
    if (maxHeterozygous > 1 || maxHeterozygous < 0) {
        stop("maxHeterozygous parameter is out of range")
    }

    # Site check
    taxa <- getTaxaList(tasObj)
    if (siteMinCount > taxa$size()) {
        stop("Minimum number of taxa exceeds total number of taxa in genotype table.")
    }

    # Range check (chromosomes)
    chroms <- getGenotypeTable(tasObj)
    chroms <- chroms$chromosomes()
    chroms <- unlist(lapply(chroms, function(x) { rJava::.jstrVal(x) }))
    if (!is.null(startChr) || !is.null(endChr)) {
        if (!(any(startChr %in% chroms)) || !(any(endChr %in% chroms))) {
            stop("Chromosome IDs not found in genotype table.")
        }
    }

    # Filter type selection
    siteRangeFilterType <- match.arg(siteRangeFilterType)
    if (missing(siteRangeFilterType) || !siteRangeFilterType %in% c("none", "sites", "position")) {
        stop(
            paste(
                "Please specify analysis type",
                "(\"none\", \"sites\", or \"position\")"
            )
        )
    }

    if (siteRangeFilterType == "position") {
        physPosLS <- getMinMaxPhysPositions(tasObj)
    }

    # Create filter siter builder plugin
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin"),
        rJava::.jnull(),
        FALSE
    )
    plugin$setParameter("siteMinCount", toString(siteMinCount))
    plugin$setParameter("siteMinAlleleFreq", toString(siteMinAlleleFreq))
    plugin$setParameter("siteMaxAlleleFreq", toString(siteMaxAlleleFreq))
    plugin$setParameter("minHeterozygous", toString(minHeterozygous))
    plugin$setParameter("maxHeterozygous", toString(maxHeterozygous))
    plugin$setParameter("removeMinorSNPStates", toString(removeMinorSNPStates))
    plugin$setParameter("removeSitesWithIndels", toString(removeSitesWithIndels))

    # Logic check necessary parameters given range filter type
    if (is.null(chrPosFile) && is.null(bedFile)) {
        if (siteRangeFilterType == "sites") {

            if (is.null(startSite) || is.null(endSite)) {
                stop("Please specify both start and end sites.")
            }

            if (endSite > jGenoTable$numberOfSites()) {
                stop("End site parameter exceeds total number of sites in genotype table.")
            }

            if (startSite > endSite) {
                stop("Start site cannot be larger than end site.")
            }

            plugin$setParameter("startSite", toString(startSite))
            plugin$setParameter("endSite", toString(endSite))

        } else if (siteRangeFilterType == "position") {

            if (is.null(startChr) || is.null(endChr)) {
                stop("Please specify both start and end chromosomes.")
            }

            if (!is.null(startPos)) startPos <- toString(startPos)
            if (!is.null(endPos)) endPos <- toString(endPos)

            if (startChr == endChr && startPos > endPos) {
                stop("Filtration paramaters outside acceptable range.")
            }

            plugin$setParameter("startChr", toString(startChr))
            plugin$setParameter("startPos", startPos)
            plugin$setParameter("endChr", toString(endChr))
            plugin$setParameter("endPos", endPos)
        }
    } else if (is.character(chrPosFile) && is.null(bedFile)) {
        tmpChrDF <- utils::read.table(chrPosFile, sep = "\t", header = TRUE)
        headCheck <- c("Chromosome", "Position")

        if (!identical(colnames(tmpChrDF), headCheck)) {
            stop("Please check chromosome position file for correct formatting")
        }

        plugin$setParameter("chrPosFile", chrPosFile)

    } else if (is.null(chrPosFile) && is.character(bedFile)) {
        plugin$setParameter("bedFile", bedFile)
    } else {
        stop("Incorrect parameter usage")
    }

    # Run plugin
    out <- tryCatch(
        {
            plugin$runPlugin(jGenoTable)
        },
        error = function(e) {
            return(-1)
        }
    )

    # Check if input had phenotype table. If yes, combine genotype with phenotype
    if (class(out) != "jobjRef") {
        message("No data returned.")
        return(NA)
    } else {
        jPhenoTable <- getPhenotypeTable(tasObj)
        if (rJava::is.jnull(jPhenoTable)) {
            .tasselObjectConstructor(out)
        } else {
            .tasselObjectConstructor(
                combineTasselGenotypePhenotype(
                    genotypeTable = out,
                    phenotype = jPhenoTable
                )
            )
        }
    }
}


#' @title Filter genotype table by taxa
#'
#' @description This function will filter R objects of
#'    \code{TasselGenotypePhenotype} class containing genotype tables.
#'    The parameters for this function are derived from TASSEL's
#'    \code{FilterTaxaBuilder} plugin.
#'
#' @name filterGenotypeTableTaxa
#' @rdname filterGenotypeTableTaxa
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#' @param minNotMissing Minimum proportion of sites not unknown to pass this
#'    filter. Value can be between 0.0 and 1.0.
#' @param minHeterozygous Minimum proportion of sites that are heterozygous.
#'    Value can be between 0.0 and 1.0.
#' @param maxHeterozygous Maximum proportion of sites that are heterozygous.
#'    Value can be between 0.0 and 1.0.
#' @param taxa Vector of taxa IDs (as character) to subset. Defaults to
#'    \code{NULL}.
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#' @export
filterGenotypeTableTaxa <- function(tasObj,
                                    minNotMissing = 0.0,
                                    minHeterozygous = 0.0,
                                    maxHeterozygous = 1.0,
                                    taxa = NULL) {

    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    # Range check
    if (minNotMissing > 1 || minNotMissing < 0 ) {
        stop("minNotMissing parameter is out of range")
    }
    if (minHeterozygous > 1 || minHeterozygous < 0 ) {
        stop("minHeterozygous parameter is out of range")
    }
    if (maxHeterozygous > 1 || maxHeterozygous < 0 ) {
        stop("maxHeterozygous parameter is out of range")
    }

    # Create filter taxa builder plugin
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin"),
        rJava::.jnull(),
        FALSE
    )
    plugin$setParameter("minNotMissing", toString(minNotMissing))
    plugin$setParameter("minHeterozygous", toString(minHeterozygous))
    plugin$setParameter("maxHeterozygous", toString(maxHeterozygous))

    if (!is.null(taxa)) {
        if (!is.vector(taxa)) {
            stop("Taxa list must be vector.")
        }
        if (!is.character(taxa)) {
            stop("Taxa list must be of type character.")
        }

        builder <- .jnew("net.maizegenetics.taxa.TaxaListBuilder")
        builder$addAll(.jarray(taxa))
        taxaArray <- builder$build()

        plugin$setParameter("includeTaxa", "true")
        plugin$setParameter("taxaList", taxaArray)
    }

    resultDataSet <- plugin$runPlugin(jGenoTable)

    # Check if input had phenotype table. If yes, combine genotype with phenotype
    jPhenoTable <- getPhenotypeTable(tasObj)
    if (rJava::is.jnull(jPhenoTable)) {
        .tasselObjectConstructor(resultDataSet)
    } else {
        .tasselObjectConstructor(
            combineTasselGenotypePhenotype(
                genotypeTable = resultDataSet,
                phenotype = jPhenoTable
            )
        )
    }
}
