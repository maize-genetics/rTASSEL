# /// Internal Helpers (verb filtering) //////////////////////////////

## ----
# Apply verb-built selectors and rebuild the class that was handed in.
#
# @param tasIn
# The list returned by '.resolveTasselInput()'.
# @param taxaSel
# A 'TaxaSelector', or NULL to leave the taxa axis alone.
# @param siteSel
# A 'SiteSelector', or NULL to leave the site axis alone.
applyVerbSelectors <- function(tasIn, taxaSel = NULL, siteSel = NULL) {
    jGt <- tasIn$jGt

    if (!is.null(taxaSel)) jGt <- applyTaxaSelector(jGt, taxaSel)
    if (!is.null(siteSel)) jGt <- applySiteSelector(jGt, siteSel)

    .wrapGenotypeResult(jGt, tasIn)
}


## ----
# Resolve a tidyselect expression against a set of IDs
#
# @param quos
# A 'quosures' list from 'rlang::enquos()'.
# @param ids
# The 'character' vector of taxa IDs or site names being selected
# from, or a 'data.frame' whose columns are the traits being selected
# from. A data frame lets predicate helpers such as 'where()' see the
# values behind each name.
# @param errorCall
# Environment of the calling verb, so that tidyselect reports failures
# against the verb rather than against this helper.
#
# @return
# An 'integer' vector of 1-based positions in 'ids'.
evalIdSelection <- function(quos, ids, errorCall) {
    data <- if (is.data.frame(ids)) {
        ids
    } else {
        stats::setNames(seq_along(ids), ids)
    }

    positions <- tidyselect::eval_select(
        rlang::expr(c(!!!quos)),
        data         = data,
        allow_rename = FALSE,
        error_call   = errorCall
    )

    unname(positions)
}


## ----
# Turn 'slice*()' positions into a keep/drop instruction
#
# Follows 'dplyr::slice()': positions are 1-based, zeros are dropped,
# positive and negative positions cannot be mixed, and positions past the
# end of the axis are ignored.
#
# @param idx
# The concatenated '...' of a slice verb.
# @param n
# Number of taxa, sites, or traits available.
# @param fn
# Name of the calling verb, used in error messages.
#
# @return
# A 'list' with an 'integer' vector 'positions' and a 'logical' 'negate',
# or NULL when nothing is left to do.
resolveSlicePositions <- function(idx, n, fn) {
    if (!is.numeric(idx)) {
        rlang::abort(c(
            sprintf("`%s()` positions must be numeric", fn),
            "x" = sprintf(
                "Got a vector of class <%s>",
                paste(class(idx), collapse = "/")
            ),
            "i" = sprintf(
                "To select by ID, use `%s()`",
                sub("^slice", "select", fn)
            )
        ))
    }

    idx <- vctrs::vec_cast(idx, integer())
    idx <- idx[!is.na(idx) & idx != 0L]

    if (any(idx > 0L) && any(idx < 0L)) {
        rlang::abort(sprintf(
            "`%s()` cannot mix positive and negative positions",
            fn
        ))
    }

    if (length(idx) == 0L) {
        return(NULL)
    }

    negate <- idx[[1L]] < 0L
    idx <- unique(abs(idx))
    idx <- idx[idx <= n]

    if (length(idx) == 0L) {
        # Every dropped position was past the end of the axis, so a
        # negative selection removes nothing at all
        if (negate) return(NULL)

        rlang::abort(c(
            sprintf("`%s()` selected nothing", fn),
            "i" = sprintf("Every position given is greater than %s", n)
        ))
    }

    list(positions = idx, negate = negate)
}


## ----
# Evaluate a verb predicate against a phenotype data mask
#
# The genotype axes hand their predicates to a selector so that simple
# thresholds can be pushed down to a TASSEL filter plugin. The phenotype
# axes have no such plugin, so they evaluate in R and this helper stands
# in for the selector: it applies the same checks and turns the mask into
# positions.
#
# @param quos
# A 'quosures' list from 'rlang::enquos()'.
# @param mask
# The data mask to evaluate against.
# @param n
# Length of the axis being filtered, which a predicate returning a
# single value is recycled to.
# @param noun
# Plural noun naming the axis, used in error messages.
#
# @return
# An 'integer' vector of 1-based positions to keep.
resolvePredicatePositions <- function(quos, mask, n, noun) {
    keep <- evalPredicate(quos, mask)

    if (!is.logical(keep)) {
        rlang::abort(sprintf(
            "A predicate over %s must evaluate to a logical vector", noun
        ))
    }

    keep <- vctrs::vec_recycle(keep, n)

    # 'NA' is not 'TRUE', so a missing value drops its row as it does in
    # 'dplyr::filter()'
    positions <- which(!is.na(keep) & keep)

    if (length(positions) == 0L) {
        rlang::abort(sprintf("No %s match the selection criteria", noun))
    }

    positions
}


## ----
# Keep whole taxa of a phenotype, however many observations they have
#
# @param tasIn
# The list returned by '.resolveTasselInput()'.
# @param ids
# A 'character' vector of taxa IDs to keep.
#
# @return
# An object of the same class as 'tasIn$original'.
keepPhenotypeTaxa <- function(tasIn, ids) {
    if (length(ids) == 0L) {
        rlang::abort("No taxa match the selection criteria")
    }

    jTaxa <- rJava::.jnew(TASSEL_JVM$TAXA_LIST_BUILDER)$
        addAll(rJava::.jarray(ids))$
        build()

    jPh <- rJava::.jnew(TASSEL_JVM$PHENO_BUILDER)$
        fromPhenotype(tasIn$jPh)$
        keepTaxa(jTaxa)$
        build()$
        get(0L)

    .wrapPhenotypeResult(jPh, tasIn)
}


## ----
# Trait attribute rows and the taxa index needed to subset them
#
# @param tasIn
# The list returned by '.resolveTasselInput()'.
#
# @return
# A 'list' with the 'tibble' 'rows' of trait attributes and the 0-based
# 'taxaIdx' of the taxa attribute.
traitAxis <- function(tasIn) {
    attrData <- phenotypeAttrData(tasIn)

    list(
        rows    = traitAttrRows(attrData),
        taxaIdx = attrData$attr_idx[attrData$trait_type == "taxa"][[1L]]
    )
}


## ----
# Apply trait positions and rebuild the class that was handed in
#
# @param tasIn
# The list returned by '.resolveTasselInput()'.
# @param axis
# The list returned by 'traitAxis()'.
# @param positions
# A 1-based 'integer' vector of positions in 'axis$rows'.
#
# @return
# An object of the same class as 'tasIn$original'.
applyTraitPositions <- function(tasIn, axis, positions) {
    if (length(positions) == 0L) {
        rlang::abort("No traits match the selection criteria")
    }

    .wrapPhenotypeResult(
        subsetPhenotypeTraits(
            tasIn$jPh,
            axis$rows$attr_idx[positions],
            axis$taxaIdx
        ),
        tasIn
    )
}



# /// Verbs (predicates) /////////////////////////////////////////////

## ----
#' @title Filter taxa, sites, or traits by metadata
#'
#' @description
#' Keeps the taxa, sites, or traits for which every given expression is
#' \code{TRUE}, in the manner of \code{dplyr::filter()}.
#' \code{filterTaxa()} and \code{filterSites()} are the pipe-friendly
#' equivalents of \code{\link{taxaWhere}()} and \code{\link{sitesWhere}()}
#' used inside \code{[}, and are backed by the same machinery, so either
#' style gives the same result.
#'
#' @details
#' On genotype data, \code{filterTaxa()} evaluates its expressions
#' against per-taxon metadata:
#'
#' \tabular{ll}{
#'    \strong{Column} \tab \strong{Description} \cr
#'    \code{taxaId} \tab Taxon name \cr
#'    \code{notMissing} \tab Proportion of sites with a genotype call \cr
#'    \code{het} \tab Proportion of sites that are heterozygous
#' }
#'
#' On phenotype data it instead evaluates them against the phenotype's
#' own columns, one element per \emph{observation}, plus a \code{taxaId}
#' alias for whichever column holds the taxa so that the same predicate
#' reads the same way on either kind of data. A phenotype may hold
#' several observations of one taxon, in which case only the matching
#' rows are kept.
#'
#' \code{filterSites()} evaluates its expressions against per-site
#' metadata:
#'
#' \tabular{ll}{
#'    \strong{Column} \tab \strong{Description} \cr
#'    \code{siteIndex} \tab 1-based position of the site in the table \cr
#'    \code{siteId} \tab Marker name \cr
#'    \code{chrom} \tab Chromosome (sequence) ID \cr
#'    \code{pos} \tab Physical position (bp) \cr
#'    \code{maf} \tab Minor allele frequency \cr
#'    \code{alleleCount} \tab Non-missing alleles (2 per diploid taxon) \cr
#'    \code{het} \tab Proportion of taxa that are heterozygous \cr
#'    \code{isIndel} \tab Whether the site contains an indel \cr
#'    \code{isBiallelic} \tab Whether the site has two or fewer states
#' }
#'
#' \code{filterTraits()} evaluates its expressions against per-trait
#' metadata:
#'
#' \tabular{ll}{
#'    \strong{Column} \tab \strong{Description} \cr
#'    \code{traitIndex} \tab 1-based position of the trait \cr
#'    \code{traitId} \tab Trait name \cr
#'    \code{traitType} \tab TASSEL attribute type: \code{"data"}, \code{"covariate"}, or \code{"factor"} \cr
#'    \code{notMissing} \tab Proportion of observations that are not missing
#' }
#'
#' The taxa column is an axis rather than a trait, so it is never
#' offered to \code{filterTraits()} and never dropped by it.
#'
#' Several expressions are combined with \code{&}, so
#' \code{filterSites(gt, maf >= 0.05, !isIndel)} and
#' \code{gt[, sitesWhere(maf >= 0.05 & !isIndel)]} are the same query. As
#' in \code{dplyr::filter()}, an expression that evaluates to \code{NA}
#' drops what it was testing. \code{\link{overlaps}()} is also available
#' inside \code{filterSites()} for range-based criteria.
#'
#' Called with no expressions, all three verbs return \code{x} unchanged.
#'
#' A \code{\linkS4class{TasselGenomicDataset}} carries both kinds of
#' data, so \code{filterTaxa()} looks at the expressions to decide which
#' to read: naming a phenotype column filters observations, and anything
#' else filters the genotype table. \code{notMissing} and \code{het}
#' therefore keep their genotype meaning throughout, and can be combined
#' with phenotype criteria in one call. Whichever axis is filtered, the
#' two components of the returned dataset are re-joined, so taxa left
#' without any observations are dropped from the genotype table as well.
#'
#' @param x An object of class \code{\linkS4class{TasselGenotype}},
#'    \code{\linkS4class{TasselPhenotype}}, or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the
#'    deprecated \code{TasselGenotypePhenotype} class are still accepted.
#'    \code{filterTraits()} needs phenotype data and
#'    \code{filterSites()} needs genotype data.
#' @param ... Expressions that evaluate to a logical vector, one element
#'    per taxon, observation, site, or trait.
#'
#' @return An object of the same class as \code{x}.
#'
#' @seealso \code{\link{taxaWhere}}, \code{\link{sitesWhere}},
#'    \code{\link{selectSites}}, \code{\link{sliceSites}},
#'    \code{\link{overlaps}}
#'
#' @examples
#' \dontrun{
#' gt <- readGenotype("path/to/genotype.hmp.txt")
#'
#' gt |> filterSites(maf >= 0.05)
#' gt |> filterSites(chrom == "1", pos >= 1e6, pos <= 2e6)
#' gt |> filterTaxa(notMissing >= 0.8, het <= 0.1)
#'
#' # Verbs chain, so both axes are filtered by piping
#' gt |>
#'     filterTaxa(startsWith(taxaId, "CML")) |>
#'     filterSites(maf >= 0.05, !isIndel)
#'
#' # On phenotype data, the trait columns are in scope
#' ph <- readPhenotype("path/to/phenotype.txt")
#'
#' ph |> filterTaxa(EarHT > 100)
#' ph |> filterTaxa(location == "A", !is.na(EarDia))
#' ph |> filterTraits(traitType == "covariate")
#' ph |> filterTraits(notMissing >= 0.95)
#'
#' # A dataset can be filtered on either kind of data, or both at once
#' ds <- readGenomicDataset(gt, ph)
#'
#' ds |> filterTaxa(notMissing >= 0.8, location == "A")
#' }
#'
#' @name filterSites
#' @rdname filterSites
#' @export
filterSites <- function(x, ...) {
    quos <- rlang::enquos(...)
    if (length(quos) == 0L) return(x)

    tasIn <- .resolveTasselInput(x, "genotype", "filterSites")

    applyVerbSelectors(tasIn, siteSel = predicateSiteSelector(quos))
}


## ----
#' @rdname filterSites
#' @export
filterTaxa <- function(x, ...) {
    quos <- rlang::enquos(...)
    if (length(quos) == 0L) return(x)

    tasIn <- .resolveTasselInput(x, "any", "filterTaxa")
    vars  <- predicateVars(quos)

    onGenotype <- function() {
        applyVerbSelectors(tasIn, taxaSel = predicateTaxaSelector(quos))
    }

    if (rJava::is.jnull(tasIn$jPh)) return(onGenotype())

    rData <- phenotypeRowData(tasIn)

    # An object carrying both kinds of data has two axes this verb could
    # filter, so the predicate decides: naming a phenotype column is a
    # phenotype query, and anything else stays on the genotype, where a
    # simple threshold can be pushed down to a TASSEL filter plugin
    onBoth <- !rJava::is.jnull(tasIn$jGt)
    if (onBoth && !any(vars %in% phenotypeOnlyVars(rData))) {
        return(onGenotype())
    }

    positions <- resolvePredicatePositions(
        quos,
        phenotypeTaxaMask(tasIn, rData, vars),
        nrow(rData),
        "observations"
    )

    .wrapPhenotypeResult(subsetPhenotypeObs(tasIn$jPh, positions), tasIn)
}


## ----
#' @rdname filterSites
#' @export
filterTraits <- function(x, ...) {
    quos <- rlang::enquos(...)
    if (length(quos) == 0L) return(x)

    tasIn <- .resolveTasselInput(x, "phenotype", "filterTraits")

    axis      <- traitAxis(tasIn)
    meta      <- buildTraitMetadata(axis$rows, phenotypeRowData(tasIn))
    positions <- resolvePredicatePositions(
        quos, meta, nrow(axis$rows), "traits"
    )

    applyTraitPositions(tasIn, axis, positions)
}



# /// Verbs (identifiers) ////////////////////////////////////////////

## ----
#' @title Select taxa, sites, or traits by ID
#'
#' @description
#' Keeps taxa, sites, or traits named by a \code{tidyselect} expression,
#' in the manner of \code{dplyr::select()}. \code{selectTaxa()} and
#' \code{selectSites()} are the pipe-friendly equivalents of
#' \code{\link{taxa}()} and \code{\link{siteIds}()} used inside \code{[},
#' with the whole \code{tidyselect} vocabulary available on top.
#'
#' @details
#' IDs can be given literally, as a vector, or with any
#' \code{tidyselect} helper, and a leading \code{-} drops rather than
#' keeps:
#'
#' \tabular{ll}{
#'    \strong{Expression} \tab \strong{Selects} \cr
#'    \code{"B73", "Mo17"} \tab Two taxa by name \cr
#'    \code{all_of(ids)} \tab Every ID in \code{ids}, erroring if any is absent \cr
#'    \code{any_of(ids)} \tab The IDs in \code{ids} that are present \cr
#'    \code{starts_with("CML")} \tab Every ID with that prefix \cr
#'    \code{matches("^PZ[AB]")} \tab Every ID matching that regular expression \cr
#'    \code{-any_of(ids)} \tab Everything except those IDs
#' }
#'
#' Bare IDs are matched exactly, so a marker name occurring more than
#' once in a table resolves to its first occurrence. Called with no
#' expressions, all three verbs return \code{x} unchanged.
#'
#' \code{selectTraits()} selects from the trait columns themselves rather
#' than from their names alone, so the \code{where()} helper can pick
#' traits out by their values, as in
#' \code{selectTraits(ph, where(is.numeric))}. The taxa column is an axis
#' rather than a trait, so it is neither selectable nor droppable.
#'
#' \code{selectTaxa()} on phenotype data works a taxon at a time rather
#' than an observation at a time: every observation of a selected taxon
#' is kept. Use \code{\link{filterTaxa}()} to select observations
#' instead.
#'
#' @param x An object of class \code{\linkS4class{TasselGenotype}},
#'    \code{\linkS4class{TasselPhenotype}}, or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the
#'    deprecated \code{TasselGenotypePhenotype} class are still accepted.
#'    \code{selectTraits()} needs phenotype data and
#'    \code{selectSites()} needs genotype data.
#' @param ... \code{tidyselect} expressions naming taxa, sites, or
#'    traits.
#'
#' @return An object of the same class as \code{x}.
#'
#' @seealso \code{\link{taxa}}, \code{\link{siteIds}},
#'    \code{\link{filterSites}}, \code{\link{sliceSites}}
#'
#' @examples
#' \dontrun{
#' gt <- readGenotype("path/to/genotype.hmp.txt")
#'
#' gt |> selectSites("PZB00859.1", "PZA01271.1")
#' gt |> selectTaxa(starts_with("CML"))
#'
#' myMarkers <- c("PZB00859.1", "PZA01271.1")
#' gt |> selectSites(all_of(myMarkers))
#'
#' # Everything but a handful of taxa
#' gt |> selectTaxa(-any_of(c("33-16", "38-11")))
#'
#' # On phenotype data, traits are selected by name or by type
#' ph <- readPhenotype("path/to/phenotype.txt")
#'
#' ph |> selectTraits(EarHT, dpoll)
#' ph |> selectTraits(starts_with("Q"))
#' ph |> selectTraits(where(is.numeric))
#' ph |> selectTraits(-EarDia)
#' }
#'
#' @name selectSites
#' @rdname selectSites
#' @export
selectSites <- function(x, ...) {
    quos <- rlang::enquos(...)
    if (length(quos) == 0L) return(x)

    tasIn <- .resolveTasselInput(x, "genotype", "selectSites")

    # Prefer the batch Java path in buildSiteMetadata() over the per-site
    # round-trips batchSiteNames() falls back to
    siteNames <- buildSiteMetadata(tasIn$jGt, needed = "siteId")$siteId
    positions <- evalIdSelection(quos, siteNames, rlang::current_env())

    if (length(positions) == 0L) {
        rlang::abort("No sites match the selection criteria")
    }

    applyVerbSelectors(tasIn, siteSel = sites(positions))
}


## ----
#' @rdname selectSites
#' @export
selectTaxa <- function(x, ...) {
    quos <- rlang::enquos(...)
    if (length(quos) == 0L) return(x)

    tasIn <- .resolveTasselInput(x, "any", "selectTaxa")

    # Selecting by ID means the same thing on either kind of data, so an
    # object carrying both is read through its genotype table
    if (rJava::is.jnull(tasIn$jGt)) {
        taxaIds   <- phenotypeTaxaNames(tasIn$jPh)
        positions <- evalIdSelection(quos, taxaIds, rlang::current_env())

        return(keepPhenotypeTaxa(tasIn, taxaIds[positions]))
    }

    taxaIds   <- batchTaxaNames(tasIn$jGt)
    positions <- evalIdSelection(quos, taxaIds, rlang::current_env())

    if (length(positions) == 0L) {
        rlang::abort("No taxa match the selection criteria")
    }

    applyVerbSelectors(tasIn, taxaSel = taxa(taxaIds[positions]))
}


## ----
#' @rdname selectSites
#' @export
selectTraits <- function(x, ...) {
    quos <- rlang::enquos(...)
    if (length(quos) == 0L) return(x)

    tasIn <- .resolveTasselInput(x, "phenotype", "selectTraits")

    axis      <- traitAxis(tasIn)
    traitData <- phenotypeRowData(tasIn)[axis$rows$trait_id]
    positions <- evalIdSelection(quos, traitData, rlang::current_env())

    applyTraitPositions(tasIn, axis, positions)
}



# /// Verbs (positions) //////////////////////////////////////////////

## ----
#' @title Select taxa, sites, or traits by position
#'
#' @description
#' Keeps taxa, sites, or traits at the given positions, in the manner of
#' \code{dplyr::slice()}. \code{sliceTaxa()} and \code{sliceSites()} are
#' the pipe-friendly equivalents of \code{\link{sites}()} and of a bare
#' numeric index used inside \code{[}.
#'
#' @details
#' Positions are 1-based, matching \code{R}'s own subsetting
#' conventions; the 0-based index TASSEL uses internally is reported in
#' the \code{Site} column of \code{\link{positionList}()}.
#'
#' As in \code{dplyr::slice()}, negative positions drop rather than keep,
#' positive and negative positions cannot be mixed, and zeros and
#' positions past the end of the axis are ignored. Called with no
#' positions, all three verbs return \code{x} unchanged.
#'
#' Taxa are counted as \code{\link{taxaList}()} reports them and traits
#' as \code{\link{traitNames}()} does, so the taxa column of a phenotype
#' is not a trait position and is never dropped.
#'
#' \code{sliceTaxa()} on phenotype data works a taxon at a time rather
#' than an observation at a time: every observation of a kept taxon is
#' kept.
#'
#' @param x An object of class \code{\linkS4class{TasselGenotype}},
#'    \code{\linkS4class{TasselPhenotype}}, or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the
#'    deprecated \code{TasselGenotypePhenotype} class are still accepted.
#'    \code{sliceTraits()} needs phenotype data and \code{sliceSites()}
#'    needs genotype data.
#' @param ... Numeric positions, all positive or all negative.
#'
#' @return An object of the same class as \code{x}.
#'
#' @seealso \code{\link{sites}}, \code{\link{selectSites}},
#'    \code{\link{filterSites}}, \code{\link{positionList}}
#'
#' @examples
#' \dontrun{
#' gt <- readGenotype("path/to/genotype.hmp.txt")
#'
#' gt |> sliceSites(1:1000)
#' gt |> sliceSites(c(10, 50, 100))
#' gt |> sliceTaxa(1:10)
#'
#' # Drop the first ten markers
#' gt |> sliceSites(-(1:10))
#'
#' # The first three traits of a phenotype, then all but the first
#' ph <- readPhenotype("path/to/phenotype.txt")
#'
#' ph |> sliceTraits(1:3)
#' ph |> sliceTraits(-1)
#' }
#'
#' @name sliceSites
#' @rdname sliceSites
#' @export
sliceSites <- function(x, ...) {
    if (...length() == 0L) return(x)

    tasIn <- .resolveTasselInput(x, "genotype", "sliceSites")

    slice <- resolveSlicePositions(
        c(...), tasIn$jGt$numberOfSites(), "sliceSites"
    )
    if (is.null(slice)) return(x)

    selector <- sites(slice$positions)
    if (slice$negate) selector <- !selector

    applyVerbSelectors(tasIn, siteSel = selector)
}


## ----
#' @rdname sliceSites
#' @export
sliceTaxa <- function(x, ...) {
    if (...length() == 0L) return(x)

    tasIn <- .resolveTasselInput(x, "any", "sliceTaxa")

    # As in 'selectTaxa()', an object carrying both kinds of data is
    # counted through its genotype table
    if (rJava::is.jnull(tasIn$jGt)) {
        taxaIds <- phenotypeTaxaNames(tasIn$jPh)

        slice <- resolveSlicePositions(c(...), length(taxaIds), "sliceTaxa")
        if (is.null(slice)) return(x)

        kept <- if (slice$negate) {
            taxaIds[-slice$positions]
        } else {
            taxaIds[slice$positions]
        }

        return(keepPhenotypeTaxa(tasIn, kept))
    }

    taxaIds <- batchTaxaNames(tasIn$jGt)

    slice <- resolveSlicePositions(c(...), length(taxaIds), "sliceTaxa")
    if (is.null(slice)) return(x)

    selector <- taxa(taxaIds[slice$positions])
    if (slice$negate) selector <- !selector

    applyVerbSelectors(tasIn, taxaSel = selector)
}


## ----
#' @rdname sliceSites
#' @export
sliceTraits <- function(x, ...) {
    if (...length() == 0L) return(x)

    tasIn <- .resolveTasselInput(x, "phenotype", "sliceTraits")

    axis  <- traitAxis(tasIn)
    slice <- resolveSlicePositions(c(...), nrow(axis$rows), "sliceTraits")
    if (is.null(slice)) return(x)

    positions <- if (slice$negate) {
        seq_len(nrow(axis$rows))[-slice$positions]
    } else {
        slice$positions
    }

    applyTraitPositions(tasIn, axis, positions)
}
