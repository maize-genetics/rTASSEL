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
# Resolve a tidyselect expression against a vector of IDs
#
# @param quos
# A 'quosures' list from 'rlang::enquos()'.
# @param ids
# The 'character' vector of taxa IDs or site names being selected from.
# @param errorCall
# Environment of the calling verb, so that tidyselect reports failures
# against the verb rather than against this helper.
#
# @return
# An 'integer' vector of 1-based positions in 'ids'.
evalIdSelection <- function(quos, ids, errorCall) {
    idPositions <- seq_along(ids)
    names(idPositions) <- ids

    positions <- tidyselect::eval_select(
        rlang::expr(c(!!!quos)),
        data         = idPositions,
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
# Number of taxa or sites available.
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
            "i" = "To select by ID, use `selectSites()` or `selectTaxa()`"
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



# /// Verbs (predicates) /////////////////////////////////////////////

## ----
#' @title Filter taxa or sites by metadata
#'
#' @description
#' Keeps the taxa or sites for which every given expression is
#' \code{TRUE}, in the manner of \code{dplyr::filter()}. These are the
#' pipe-friendly equivalents of \code{\link{taxaWhere}()} and
#' \code{\link{sitesWhere}()} used inside \code{[}, and are backed by the
#' same machinery, so either style gives the same result.
#'
#' @details
#' \code{filterTaxa()} evaluates its expressions against per-taxon
#' metadata:
#'
#' \tabular{ll}{
#'    \strong{Column} \tab \strong{Description} \cr
#'    \code{taxaId} \tab Taxon name \cr
#'    \code{notMissing} \tab Proportion of sites with a genotype call \cr
#'    \code{het} \tab Proportion of sites that are heterozygous
#' }
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
#' Several expressions are combined with \code{&}, so
#' \code{filterSites(gt, maf >= 0.05, !isIndel)} and
#' \code{gt[, sitesWhere(maf >= 0.05 & !isIndel)]} are the same query.
#' \code{\link{overlaps}()} is also available inside
#' \code{filterSites()} for range-based criteria.
#'
#' Called with no expressions, both verbs return \code{x} unchanged.
#'
#' @param x An object of class \code{\linkS4class{TasselGenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the
#'    deprecated \code{TasselGenotypePhenotype} class are still accepted.
#' @param ... Expressions that evaluate to a logical vector, one element
#'    per taxon or site.
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

    tasIn <- .resolveTasselInput(x, "genotype", "filterTaxa")

    applyVerbSelectors(tasIn, taxaSel = predicateTaxaSelector(quos))
}



# /// Verbs (identifiers) ////////////////////////////////////////////

## ----
#' @title Select taxa or sites by ID
#'
#' @description
#' Keeps taxa or sites named by a \code{tidyselect} expression, in the
#' manner of \code{dplyr::select()}. These are the pipe-friendly
#' equivalents of \code{\link{taxa}()} and \code{\link{siteIds}()} used
#' inside \code{[}, with the whole \code{tidyselect} vocabulary available
#' on top.
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
#' expressions, both verbs return \code{x} unchanged.
#'
#' @param x An object of class \code{\linkS4class{TasselGenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the
#'    deprecated \code{TasselGenotypePhenotype} class are still accepted.
#' @param ... \code{tidyselect} expressions naming taxa or sites.
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

    tasIn <- .resolveTasselInput(x, "genotype", "selectTaxa")

    taxaIds   <- batchTaxaNames(tasIn$jGt)
    positions <- evalIdSelection(quos, taxaIds, rlang::current_env())

    if (length(positions) == 0L) {
        rlang::abort("No taxa match the selection criteria")
    }

    applyVerbSelectors(tasIn, taxaSel = taxa(taxaIds[positions]))
}



# /// Verbs (positions) //////////////////////////////////////////////

## ----
#' @title Select taxa or sites by position
#'
#' @description
#' Keeps taxa or sites at the given positions, in the manner of
#' \code{dplyr::slice()}. These are the pipe-friendly equivalents of
#' \code{\link{sites}()} and of a bare numeric index used inside
#' \code{[}.
#'
#' @details
#' Positions are 1-based, matching \code{R}'s own subsetting
#' conventions; the 0-based index TASSEL uses internally is reported in
#' the \code{Site} column of \code{\link{positionList}()}.
#'
#' As in \code{dplyr::slice()}, negative positions drop rather than keep,
#' positive and negative positions cannot be mixed, and zeros and
#' positions past the end of the axis are ignored. Called with no
#' positions, both verbs return \code{x} unchanged.
#'
#' @param x An object of class \code{\linkS4class{TasselGenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the
#'    deprecated \code{TasselGenotypePhenotype} class are still accepted.
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

    tasIn <- .resolveTasselInput(x, "genotype", "sliceTaxa")

    taxaIds <- batchTaxaNames(tasIn$jGt)

    slice <- resolveSlicePositions(c(...), length(taxaIds), "sliceTaxa")
    if (is.null(slice)) return(x)

    selector <- taxa(taxaIds[slice$positions])
    if (slice$negate) selector <- !selector

    applyVerbSelectors(tasIn, taxaSel = selector)
}
