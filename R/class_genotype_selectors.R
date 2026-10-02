# /// Selector Classes ///////////////////////////////////////////////

## ----
#' @title TaxaSelector Class
#'
#' @description
#' S4 class representing taxa selection criteria for bracket-based
#' filtering of \code{TasselGenotype}, \code{TasselPhenotype}, and
#' \code{TasselGenomicDataset} objects.
#'
#' @slot type Character indicating selector type: \code{"ids"} or
#'   \code{"predicate"}.
#' @slot ids Character vector of taxa IDs (used when
#'   \code{type = "ids"}).
#' @slot quo A quosure, or a list of quosures combined with \code{&},
#'   for predicate evaluation (used when \code{type = "predicate"}).
#' @slot negate Logical indicating whether to negate the selection.
#'
#' @name TaxaSelector-class
#' @rdname TaxaSelector-class
#' @exportClass TaxaSelector
setClass("TaxaSelector", slots = c(
    type   = "character",
    ids    = "character",
    quo    = "ANY",
    negate = "logical"
))

## ----
#' @title SiteSelector Class
#'
#' @description
#' S4 class representing site selection criteria for bracket-based
#' filtering of \code{TasselGenotype} and \code{TasselGenomicDataset}
#' objects.
#'
#' @slot type Character indicating selector type: \code{"indices"},
#'   \code{"names"}, \code{"chrom"}, \code{"region"},
#'   \code{"granges"}, or \code{"predicate"}.
#' @slot indices Integer vector of 1-based site indices (used when
#'   \code{type = "indices"}).
#' @slot ids Character vector of site/marker names (used when
#'   \code{type = "names"}).
#' @slot chromId Character vector of chromosome IDs (used when
#'   \code{type = "chrom"} or \code{"region"}).
#' @slot start Numeric start position in bp (used when
#'   \code{type = "region"}).
#' @slot end Numeric end position in bp (used when
#'   \code{type = "region"}).
#' @slot granges A \code{GRanges} object (used when
#'   \code{type = "granges"}).
#' @slot quo A quosure, or a list of quosures combined with \code{&},
#'   for predicate evaluation (used when \code{type = "predicate"}).
#' @slot negate Logical indicating whether to negate the selection.
#'
#' @name SiteSelector-class
#' @rdname SiteSelector-class
#' @exportClass SiteSelector
setClass("SiteSelector", slots = c(
    type    = "character",
    indices = "integer",
    ids     = "character",
    chromId = "character",
    start   = "numeric",
    end     = "numeric",
    granges = "ANY",
    quo     = "ANY",
    negate  = "logical"
))


# /// Selector Constructors //////////////////////////////////////////

## ----
#' @title Select Taxa by ID
#'
#' @description
#' Creates a \code{\linkS4class{TaxaSelector}} for filtering by taxa
#' names.
#'
#' @details
#' On phenotype data a taxon may hold more than one observation, in
#' which case every observation of a named taxon is kept.
#'
#' @param ... Character taxa IDs to select.
#'
#' @return A \code{\linkS4class{TaxaSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[taxa("B73", "Mo17"), ]
#' ph[taxa("B73", "Mo17"), ]
#' }
#'
#' @export
taxa <- function(...) {
    ids <- as.character(c(...))
    if (length(ids) == 0) rlang::abort("At least one taxon ID must be provided")
    methods::new("TaxaSelector",
        type = "ids", ids = ids, quo = NULL, negate = FALSE
    )
}

## ----
#' @title Select Taxa by Predicate
#'
#' @description
#' Creates a \code{\linkS4class{TaxaSelector}} using a predicate
#' expression evaluated against taxa metadata. On a genotype table
#' the available columns in the data mask are \code{taxaId} (taxon
#' name), \code{notMissing} (proportion of sites with a genotype
#' call), and \code{het} (proportion of sites that are
#' heterozygous).
#'
#' @details
#' On phenotype data the mask is the phenotype's own columns, one
#' element per observation, plus a \code{taxaId} alias for whichever
#' column holds the taxa, so the same predicate reads the same way on
#' either kind of data. \code{notMissing} and \code{het} are genotype
#' metrics and are not in scope on a phenotype by itself.
#'
#' A \code{\linkS4class{TasselGenomicDataset}} carries both, so the
#' predicate picks the axis: naming a phenotype column filters
#' observations, and anything else filters the taxa of the genotype
#' table. This is the same rule \code{\link{filterTaxa}()} applies.
#'
#' @param expr An unquoted expression evaluated against taxa metadata.
#'
#' @return A \code{\linkS4class{TaxaSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[taxaWhere(startsWith(taxaId, "NAM")), ]
#' gt[taxaWhere(notMissing >= 0.8), ]
#' gt[taxaWhere(het <= 0.1), ]
#'
#' # On phenotype data the trait columns are in scope
#' ph[taxaWhere(EarHT > 100), ]
#' ph[taxaWhere(location == "A" & !is.na(EarDia)), ]
#' }
#'
#' @export
taxaWhere <- function(expr) {
    predicateTaxaSelector(rlang::enquo(expr))
}

## ----
# Wrap a predicate in a TaxaSelector.  Shared by taxaWhere(), which
# passes a single quosure, and filterTaxa(), which passes a list of them.
predicateTaxaSelector <- function(quo) {
    methods::new("TaxaSelector",
        type = "predicate", ids = character(0), quo = quo, negate = FALSE
    )
}

## ----
#' @title Select Sites by Index
#'
#' @description
#' Creates a \code{\linkS4class{SiteSelector}} for filtering by site
#' index. Indices are 1-based, matching R's own subsetting
#' conventions; the 0-based index TASSEL uses internally is reported
#' in the \code{Site} column of \code{\link{positionList}()}.
#'
#' @param ... Integer site indices (1-based).
#'
#' @return A \code{\linkS4class{SiteSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[, sites(1:1000)]
#' gt[, sites(c(10, 50, 100))]
#' }
#'
#' @export
sites <- function(...) {
    idx <- as.integer(c(...))
    if (length(idx) == 0) rlang::abort("At least one site index must be provided")
    if (any(idx < 1)) {
        rlang::abort(c(
            "Site indices must be 1-based",
            "x" = "Got an index less than 1"
        ))
    }
    methods::new("SiteSelector",
        type = "indices", indices = idx, ids = character(0),
        chromId = character(0), start = numeric(0), end = numeric(0),
        granges = NULL, quo = NULL, negate = FALSE
    )
}

## ----
#' @title Select Sites by Name
#'
#' @description
#' Creates a \code{\linkS4class{SiteSelector}} for filtering by
#' marker or SNP ID strings.
#'
#' @param ... Character site name strings.
#'
#' @return A \code{\linkS4class{SiteSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[, siteIds("rs1", "rs2")]
#' }
#'
#' @export
siteIds <- function(...) {
    ids <- as.character(c(...))
    if (length(ids) == 0) rlang::abort("At least one site ID must be provided")
    methods::new("SiteSelector",
        type = "names", indices = integer(0), ids = ids,
        chromId = character(0), start = numeric(0), end = numeric(0),
        granges = NULL, quo = NULL, negate = FALSE
    )
}

## ----
#' @title Select Sites by Chromosome
#'
#' @description
#' Creates a \code{\linkS4class{SiteSelector}} for filtering by one
#' or more chromosome IDs.
#'
#' @param ... Character chromosome IDs.
#'
#' @return A \code{\linkS4class{SiteSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[, chrom("chr3")]
#' gt[, chrom("1", "5", "10")]
#' }
#'
#' @export
chrom <- function(...) {
    chromIds <- as.character(c(...))
    if (length(chromIds) == 0) rlang::abort("At least one chromosome ID must be provided")
    methods::new("SiteSelector",
        type = "chrom", indices = integer(0), ids = character(0),
        chromId = chromIds, start = numeric(0), end = numeric(0),
        granges = NULL, quo = NULL, negate = FALSE
    )
}

## ----
#' @title Select Sites by Genomic Region
#'
#' @description
#' Creates a \code{\linkS4class{SiteSelector}} for filtering by
#' genomic coordinates, given either as a single chromosome and
#' coordinate range or as a \code{GRanges} object holding any number
#' of ranges.
#'
#' @details
#' A \code{GRanges} object is the migration path for the
#' \code{bedFile} and \code{chrPosFile} arguments of the deprecated
#' \code{\link{filterGenotypeTableSites}()}: read a BED file with
#' \code{rtracklayer::import()} and build ranges from a chromosome
#' and position table with
#' \code{GenomicRanges::GRanges(chrom, IRanges::IRanges(pos, pos))}.
#'
#' @param x A character chromosome ID, or a \code{GRanges} object. If
#'   a \code{GRanges} is supplied, \code{start} and \code{end} must be
#'   omitted.
#' @param start Numeric start position in bp.
#' @param end Numeric end position in bp.
#'
#' @return A \code{\linkS4class{SiteSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[, region("chr1", 1e6, 2e6)]
#'
#' gr <- GenomicRanges::GRanges(
#'     seqnames = c("chr1", "chr2"),
#'     ranges   = IRanges::IRanges(start = c(1e6, 5e5), end = c(2e6, 1e6))
#' )
#' gt[, region(gr)]
#' }
#'
#' @export
region <- function(x, start, end) {
    if (methods::is(x, "GRanges")) {
        if (!missing(start) || !missing(end)) {
            rlang::abort(c(
                "`start` and `end` cannot be used with a <GRanges> object",
                "i" = "Encode the coordinates in the <GRanges> object instead"
            ))
        }
        if (length(x) == 0) {
            rlang::abort("At least one range must be provided")
        }

        return(methods::new("SiteSelector",
            type = "granges", indices = integer(0), ids = character(0),
            chromId = character(0), start = numeric(0), end = numeric(0),
            granges = x, quo = NULL, negate = FALSE
        ))
    }

    methods::new("SiteSelector",
        type = "region", indices = integer(0), ids = character(0),
        chromId = as.character(x), start = as.numeric(start),
        end = as.numeric(end), granges = NULL, quo = NULL, negate = FALSE
    )
}

## ----
#' @title Select Sites by Predicate
#'
#' @description
#' Creates a \code{\linkS4class{SiteSelector}} using a predicate
#' expression evaluated against site metadata. Available columns in
#' the data mask: \code{siteIndex} (1-based), \code{siteId},
#' \code{chrom}, \code{pos}, \code{maf}, \code{alleleCount},
#' \code{het}, \code{isIndel}, \code{isBiallelic}.
#'
#' @param expr An unquoted expression evaluated against site metadata.
#'
#' @return A \code{\linkS4class{SiteSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[, sitesWhere(maf >= 0.05)]
#' gt[, sitesWhere(chrom == "chr1" & maf >= 0.05)]
#' gt[, sitesWhere(alleleCount >= 10)]
#' gt[, sitesWhere(het <= 0.5)]
#' gt[, sitesWhere(!isIndel)]
#' gt[, sitesWhere(isBiallelic)]
#' }
#'
#' @export
sitesWhere <- function(expr) {
    predicateSiteSelector(rlang::enquo(expr))
}

## ----
# Wrap a predicate in a SiteSelector.  Shared by sitesWhere(), which
# passes a single quosure, and filterSites(), which passes a list of them.
predicateSiteSelector <- function(quo) {
    methods::new("SiteSelector",
        type = "predicate", indices = integer(0), ids = character(0),
        chromId = character(0), start = numeric(0), end = numeric(0),
        granges = NULL, quo = quo, negate = FALSE
    )
}


## ----
#' @title Select Sites Overlapping Genomic Ranges
#'
#' @description
#' Tests whether each site falls inside a set of genomic ranges. This
#' function has no use on its own: it is only meaningful inside
#' \code{\link{sitesWhere}()} or \code{\link{filterSites}()}, where the
#' site metadata it needs is in scope. It is the predicate-friendly
#' counterpart of \code{\link{region}()}, and unlike \code{region()} it
#' can be combined with other site criteria in one expression.
#'
#' @param ranges A \code{GRanges} object.
#'
#' @return A \code{logical} vector with one element per site.
#'
#' @seealso \code{\link{region}}, \code{\link{sitesWhere}},
#'    \code{\link{filterSites}}
#'
#' @examples
#' \dontrun{
#' gr <- GenomicRanges::GRanges(
#'     seqnames = c("1", "2"),
#'     ranges   = IRanges::IRanges(start = c(1e6, 5e5), end = c(2e6, 1e6))
#' )
#'
#' gt[, sitesWhere(overlaps(gr))]
#' gt |> filterSites(overlaps(gr), maf >= 0.05)
#' }
#'
#' @export
overlaps <- function(ranges) {
    rlang::abort(c(
        "`overlaps()` must be used inside `sitesWhere()` or `filterSites()`",
        "i" = "Outside a site predicate there is no site metadata to test"
    ))
}


# /// Negation Methods ///////////////////////////////////////////////

## ----
#' @param x A \code{TaxaSelector} object.
#'
#' @rdname TaxaSelector-class
#' @aliases !,TaxaSelector-method
setMethod("!", "TaxaSelector", function(x) {
    x@negate <- !x@negate
    x
})

## ----
#' @param x A \code{SiteSelector} object.
#'
#' @rdname SiteSelector-class
#' @aliases !,SiteSelector-method
setMethod("!", "SiteSelector", function(x) {
    x@negate <- !x@negate
    x
})



