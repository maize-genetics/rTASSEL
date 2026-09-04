## ----
#' @title Filter genotype table by sites
#'
#' @description
#' \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
#'
#' This function will filter R objects containing genotype tables by
#' marker (site) criteria. The parameters for this function are derived
#' from TASSEL's \code{FilterSiteBuilder} plugin.
#'
#' Use bracket subsetting with the site selectors instead
#' (\code{gt[, <selector>]}); see \emph{Details} for the equivalent of
#' each argument.
#'
#' @details
#' Every argument of this function has a bracket-API replacement:
#'
#' \tabular{ll}{
#'    \strong{Deprecated argument} \tab \strong{Replacement} \cr
#'    \code{siteMinCount = 150} \tab \code{gt[, sitesWhere(alleleCount >= 150)]} \cr
#'    \code{siteMinAlleleFreq = 0.05} \tab \code{gt[, sitesWhere(maf >= 0.05)]} \cr
#'    \code{siteMaxAlleleFreq = 0.3} \tab \code{gt[, sitesWhere(maf <= 0.3)]} \cr
#'    \code{minHeterozygous = 0.1} \tab \code{gt[, sitesWhere(het >= 0.1)]} \cr
#'    \code{maxHeterozygous = 0.1} \tab \code{gt[, sitesWhere(het <= 0.1)]} \cr
#'    \code{removeSitesWithIndels = TRUE} \tab \code{gt[, sitesWhere(!isIndel)]} \cr
#'    \code{removeMinorSNPStates = TRUE} \tab \code{removeMinorSNPStates(gt)} \cr
#'    \code{siteRangeFilterType = "sites"} \tab \code{gt[, sites(startSite:endSite)]} \cr
#'    \code{siteRangeFilterType = "position"} \tab \code{gt[, region("1", startPos, endPos)]} \cr
#'    \code{gRangesObj = gr} \tab \code{gt[, region(gr)]} \cr
#'    \code{bedFile = "x.bed"} \tab \code{gt[, region(rtracklayer::import("x.bed"))]} \cr
#'    \code{chrPosFile = "x.tsv"} \tab \code{gt[, region(gr)]}, building \code{gr} from the file
#' }
#'
#' Site indices passed to \code{\link{sites}()} are 1-based, while
#' \code{startSite} and \code{endSite} are the 0-based indices TASSEL uses
#' internally: \code{startSite = 1, endSite = 3} becomes
#' \code{gt[, sites(2:4)]}.
#'
#' Unlike this function, selectors can be combined with \code{&} inside a
#' single \code{\link{sitesWhere}()} call and taxa can be filtered in the
#' same expression:
#' \code{gt[taxaWhere(notMissing >= 0.8), sitesWhere(maf >= 0.05 & !isIndel)]}.
#'
#' @name filterGenotypeTableSites
#' @rdname filterGenotypeTableSites
#'
#' @param tasObj An object of class \code{\linkS4class{TasselGenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the deprecated
#'    \code{TasselGenotypePhenotype} class are still accepted.
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
#' @param gRangesObj Filter genotype table by \code{GenomicRanges} object.
#'    If this parameter is selected, you cannot utilize the parameters,
#'    \code{chrPosFile} or \code{bedFile}. Defaults to \code{NULL}.
#' @param chrPosFile An optional chromosome position file path of
#'    \code{character} class. Defaults to \code{NULL}. \strong{Note:}
#'    a chromosome position file must contain correct formatting
#'    (e.g. a two column file with the header of
#'    \code{c("Chromosome", "Position")}).
#' @param bedFile An optional BED coordinate file path of
#'    \code{character} class. Defaults to \code{NULL}.
#'
#' @return Returns an object of the same class as \code{tasObj}.
#'
#' @seealso \code{\link{sitesWhere}}, \code{\link{sites}},
#'    \code{\link{siteIds}}, \code{\link{chrom}}, \code{\link{region}},
#'    \code{\link{removeMinorSNPStates}},
#'    \code{\linkS4class{TasselGenotype}}
#'
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom rJava .jarray
#' @importFrom rJava .jnull
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @export
filterGenotypeTableSites <- function(
    tasObj,
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
    gRangesObj = NULL,
    chrPosFile = NULL,
    bedFile = NULL
) {
    lifecycle::deprecate_warn(
        when    = "0.14.0",
        what    = "filterGenotypeTableSites()",
        details = c(
            "i" = "Use bracket subsetting instead, e.g. `gt[, sitesWhere(maf >= 0.05)]`.",
            "i" = "See `?filterGenotypeTableSites` for a full argument mapping."
        )
    )

    # The notice above already flags legacy usage, so the one
    # '.resolveTasselInput()' raises for <TasselGenotypePhenotype>
    # input would only repeat it
    rlang::local_options(lifecycle_verbosity = "quiet")

    tasIn <- .resolveTasselInput(
        tasObj, "genotype", "filterGenotypeTableSites"
    )
    jGenoTable <- tasIn$jGt

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
    jTaxa <- getTaxaList(tasObj)
    if (siteMinCount > jTaxa$size()) {
        stop("Minimum number of taxa exceeds total number of taxa in genotype table.")
    }

    # Range check (chromosomes)
    chroms <- jGenoTable$chromosomes()
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

            if (startSite < 0 || endSite < 0) {
                stop("startSite and endSite must be non-negative integers.")
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

            if (!is.null(startPos) && startPos < 0) {
                stop("startPos must be a non-negative integer.")
            }
            if (!is.null(endPos) && endPos < 0) {
                stop("endPos must be a non-negative integer.")
            }

            if (startChr == endChr && !is.null(startPos) && !is.null(endPos) && startPos > endPos) {
                stop("Filtration parameters outside acceptable range.")
            }

            if (!is.null(startPos)) startPos <- as.character(startPos)
            if (!is.null(endPos)) endPos <- as.character(endPos)

            plugin$setParameter("startChr", toString(startChr))
            plugin$setParameter("startPos", startPos)
            plugin$setParameter("endChr", toString(endChr))
            plugin$setParameter("endPos", endPos)
        }
    } else if (is.character(chrPosFile) && is.null(bedFile) && is.null(gRangesObj)) {
        tmpChrDF <- utils::read.table(chrPosFile, sep = "\t", header = TRUE)
        headCheck <- c("Chromosome", "Position")

        if (!identical(colnames(tmpChrDF), headCheck)) {
            stop("Please check chromosome position file for correct formatting")
        }

        plugin$setParameter("chrPosFile", chrPosFile)

    } else if (is.null(chrPosFile) && is.character(bedFile) && is.null(gRangesObj)) {
        plugin$setParameter("bedFile", bedFile)
    } else {
        stop("Incorrect parameter usage")
    }

    # Filter by GRanges object if present
    if (!is.null(gRangesObj)) {
        if (!inherits(gRangesObj, "GRanges")) {
            stop("gRangesObj must be of class `GRanges`")
        }
        jRC <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
        jGenoTable <- jRC$filterSitesByGRanges(
            jGenoTable,
            rJava::.jarray(as.vector(GenomicRanges::seqnames(gRangesObj))),
            rJava::.jarray(GenomicRanges::start(gRangesObj)),
            rJava::.jarray(GenomicRanges::end(gRangesObj))
        )
    }

    # Try to run plugin
    out <- tryCatch(
        {
            plugin$runPlugin(jGenoTable)
        },
        error = function(e) {
            return(-1)
        }
    )

    if (!inherits(out, "jobjRef")) {
        message("No data returned.")
        return(NA)
    }

    .wrapGenotypeResult(out, tasIn)
}


## ----
#' @title Filter genotype table by site IDs
#'
#' @description
#' \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
#'
#' Filter a genotype table object by specifying literal site names (IDs)
#' for variant markers.
#'
#' Use bracket subsetting with \code{\link{siteIds}()} instead; see
#' \emph{Details}.
#'
#' @details
#' \tabular{ll}{
#'    \strong{Deprecated argument} \tab \strong{Replacement} \cr
#'    \code{siteNames = c("rs1", "rs2")} \tab \code{gt[, siteIds("rs1", "rs2")]}
#' }
#'
#' A character vector can be handed to \code{\link{siteIds}()} directly
#' (\code{gt[, siteIds(myMarkers)]}), or passed as the \code{j} index
#' without a selector (\code{gt[, myMarkers]}).
#'
#' @name filterGenotypeTableBySiteName
#' @rdname filterGenotypeTableBySiteName
#'
#' @param tasObj An object of class \code{\linkS4class{TasselGenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the deprecated
#'    \code{TasselGenotypePhenotype} class are still accepted.
#' @param siteNames A character vector of site names to filter on.
#'
#' @return Returns an object of the same class as \code{tasObj}.
#'
#' @seealso \code{\link{siteIds}}, \code{\linkS4class{TasselGenotype}}
#'
#' @export
filterGenotypeTableBySiteName <- function(tasObj, siteNames) {
    lifecycle::deprecate_warn(
        when    = "0.14.0",
        what    = "filterGenotypeTableBySiteName()",
        details = c(
            "i" = 'Use bracket subsetting instead, e.g. `gt[, siteIds("rs1", "rs2")]`.'
        )
    )

    rlang::local_options(lifecycle_verbosity = "quiet")

    tasIn <- .resolveTasselInput(
        tasObj, "genotype", "filterGenotypeTableBySiteName"
    )
    gtJava <- tasIn$jGt

    # Instantiate new ArrayList object and populate with site names
    idsToKeep <- rJava::.jnew("java.util.ArrayList")
    for (id in siteNames) idsToKeep$add(id)

    # Instantiate filter site builder plugin
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin"),
        rJava::.jnull(),
        FALSE
    )
    plugin$siteNamesList(idsToKeep)

    # Instantiate DataSet object
    dataSet <- rJava::J("net.maizegenetics.plugindef.DataSet")

    # Try to filter data - if outofbounds exception - no sites returned
    gtFilter <- tryCatch(
        expr = {
            .wrapGenotypeResult(
                plugin$runPlugin(dataSet$getDataSet(gtJava)),
                tasIn
            )
        },
        error = function(e) {
            message("Error. No sites found: ")
            message(" ", conditionMessage(e), "\n")
        }
    )

    return(gtFilter)
}


## ----
#' @title Filter genotype table by taxa
#'
#' @description
#' \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
#'
#' This function will filter R objects containing genotype tables by taxa
#' (sample) criteria. The parameters for this function are derived from
#' TASSEL's \code{FilterTaxaBuilder} plugin.
#'
#' Use bracket subsetting with the taxa selectors instead
#' (\code{gt[<selector>, ]}); see \emph{Details} for the equivalent of each
#' argument.
#'
#' @details
#' Every argument of this function has a bracket-API replacement:
#'
#' \tabular{ll}{
#'    \strong{Deprecated argument} \tab \strong{Replacement} \cr
#'    \code{minNotMissing = 0.8} \tab \code{gt[taxaWhere(notMissing >= 0.8), ]} \cr
#'    \code{minHeterozygous = 0.1} \tab \code{gt[taxaWhere(het >= 0.1), ]} \cr
#'    \code{maxHeterozygous = 0.1} \tab \code{gt[taxaWhere(het <= 0.1), ]} \cr
#'    \code{taxa = c("B73", "Mo17")} \tab \code{gt[taxa("B73", "Mo17"), ]}
#'  }
#'
#' \code{\link{taxaWhere}()} also evaluates arbitrary R expressions against
#' the \code{taxaId} column, which has no equivalent here:
#' \code{gt[taxaWhere(startsWith(taxaId, "NAM")), ]}.
#'
#' @name filterGenotypeTableTaxa
#' @rdname filterGenotypeTableTaxa
#'
#' @param tasObj An object of class \code{\linkS4class{TasselGenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the deprecated
#'    \code{TasselGenotypePhenotype} class are still accepted.
#' @param minNotMissing Minimum proportion of sites not unknown to pass this
#'    filter. Value can be between 0.0 and 1.0.
#' @param minHeterozygous Minimum proportion of sites that are heterozygous.
#'    Value can be between 0.0 and 1.0.
#' @param maxHeterozygous Maximum proportion of sites that are heterozygous.
#'    Value can be between 0.0 and 1.0.
#' @param taxa Vector of taxa IDs (as character) to subset. Defaults to
#'    \code{NULL}.
#'
#' @return Returns an object of the same class as \code{tasObj}.
#'
#' @seealso \code{\link{taxaWhere}}, \code{\link{taxa}},
#'    \code{\linkS4class{TasselGenotype}}
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#' @export
filterGenotypeTableTaxa <- function(
    tasObj,
    minNotMissing = 0.0,
    minHeterozygous = 0.0,
    maxHeterozygous = 1.0,
    taxa = NULL
) {
    lifecycle::deprecate_warn(
        when    = "0.14.0",
        what    = "filterGenotypeTableTaxa()",
        details = c(
            "i" = "Use bracket subsetting instead, e.g. `gt[taxaWhere(notMissing >= 0.8), ]`.",
            "i" = "See `?filterGenotypeTableTaxa` for a full argument mapping."
        )
    )

    rlang::local_options(lifecycle_verbosity = "quiet")

    tasIn <- .resolveTasselInput(
        tasObj, "genotype", "filterGenotypeTableTaxa"
    )
    jGenoTable <- tasIn$jGt

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

    .wrapGenotypeResult(plugin$runPlugin(jGenoTable), tasIn)
}
