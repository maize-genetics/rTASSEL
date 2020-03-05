#---------------------------------------------------------------------
# Script Name:   FilterFunctions.R
# Description:   Functions for TASSEL filtration
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2019-04-04 at 21:15:54
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
#' @param siteRangeFilterType True if filtering by site numbers. False if
#'    filtering by chromosome and position. Options are
#'    \code{NONE}, \code{SITES}, or \code{POSITIONS}. Defaults to \code{NONE}.
#' @param startSite The start site. Defaults to 0.
#' @param endSite The end site. Defaults to 0.
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
                                     siteRangeFilterType = "NONE",
                                     startSite = 0,
                                     endSite = 0) {

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
    plugin$setParameter("siteRangeFilterType", toString(siteRangeFilterType))
    plugin$setParameter("startSite", toString(startSite))
    plugin$setParameter("endSite", toString(endSite))
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
                                    maxHeterozygous = 1.0) {

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
    # plugin$setParameter("includeTaxa", toString(includeTaxa))
    # plugin$taxaList(taxaList)
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
