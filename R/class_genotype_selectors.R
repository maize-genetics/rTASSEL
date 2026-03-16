# /// Selector Classes ///////////////////////////////////////////////

## ----
#' @title TaxaSelector Class
#'
#' @description
#' S4 class representing taxa selection criteria for bracket-based
#' filtering of \code{TasselGenotype} objects.
#'
#' @slot type Character indicating selector type: \code{"ids"} or
#'   \code{"predicate"}.
#' @slot ids Character vector of taxa IDs (used when
#'   \code{type = "ids"}).
#' @slot quo Quosure for predicate evaluation (used when
#'   \code{type = "predicate"}).
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
#' filtering of \code{TasselGenotype} objects.
#'
#' @slot type Character indicating selector type: \code{"indices"},
#'   \code{"names"}, \code{"chrom"}, \code{"region"}, or
#'   \code{"predicate"}.
#' @slot indices Integer vector of 0-based TASSEL site indices (used
#'   when \code{type = "indices"}).
#' @slot ids Character vector of site/marker names (used when
#'   \code{type = "names"}).
#' @slot chromId Character vector of chromosome IDs (used when
#'   \code{type = "chrom"} or \code{"region"}).
#' @slot start Numeric start position in bp (used when
#'   \code{type = "region"}).
#' @slot end Numeric end position in bp (used when
#'   \code{type = "region"}).
#' @slot quo Quosure for predicate evaluation (used when
#'   \code{type = "predicate"}).
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
    quo     = "ANY",
    negate  = "logical"
))


# /// Selector Constructors //////////////////////////////////////////

## ----
#' @title Select Taxa by ID
#'
#' @description
#' Creates a \code{\linkS4class{TaxaSelector}} for filtering a
#' \code{TasselGenotype} by taxa names.
#'
#' @param ... Character taxa IDs to select.
#'
#' @return A \code{\linkS4class{TaxaSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[taxa("B73", "Mo17"), ]
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
#' expression evaluated against taxa metadata.
#'
#' @param expr An unquoted expression evaluated against a data frame
#'   containing a \code{taxaId} column.
#'
#' @return A \code{\linkS4class{TaxaSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[taxaWhere(startsWith(taxaId, "NAM")), ]
#' }
#'
#' @export
taxaWhere <- function(expr) {
    quo <- rlang::enquo(expr)
    methods::new("TaxaSelector",
        type = "predicate", ids = character(0), quo = quo, negate = FALSE
    )
}

## ----
#' @title Select Sites by Index
#'
#' @description
#' Creates a \code{\linkS4class{SiteSelector}} for filtering by
#' 0-based TASSEL site indices.
#'
#' @param ... Integer site indices (0-based).
#'
#' @return A \code{\linkS4class{SiteSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[, sites(0:999)]
#' gt[, sites(c(10, 50, 100))]
#' }
#'
#' @export
sites <- function(...) {
    idx <- as.integer(c(...))
    if (length(idx) == 0) rlang::abort("At least one site index must be provided")
    methods::new("SiteSelector",
        type = "indices", indices = idx, ids = character(0),
        chromId = character(0), start = numeric(0), end = numeric(0),
        quo = NULL, negate = FALSE
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
        quo = NULL, negate = FALSE
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
        quo = NULL, negate = FALSE
    )
}

## ----
#' @title Select Sites by Genomic Region
#'
#' @description
#' Creates a \code{\linkS4class{SiteSelector}} for filtering by a
#' chromosomal coordinate range.
#'
#' @param chrom Character chromosome ID.
#' @param start Numeric start position in bp.
#' @param end Numeric end position in bp.
#'
#' @return A \code{\linkS4class{SiteSelector}} object.
#'
#' @examples
#' \dontrun{
#' gt[, region("chr1", 1e6, 2e6)]
#' }
#'
#' @export
region <- function(chrom, start, end) {
    methods::new("SiteSelector",
        type = "region", indices = integer(0), ids = character(0),
        chromId = as.character(chrom), start = as.numeric(start),
        end = as.numeric(end), quo = NULL, negate = FALSE
    )
}

## ----
#' @title Select Sites by Predicate
#'
#' @description
#' Creates a \code{\linkS4class{SiteSelector}} using a predicate
#' expression evaluated against site metadata. Available columns in
#' the data mask: \code{siteIndex}, \code{siteId}, \code{chrom},
#' \code{pos}, \code{maf}, \code{alleleCount}, \code{het},
#' \code{isIndel}, \code{isBiallelic}.
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
    quo <- rlang::enquo(expr)
    methods::new("SiteSelector",
        type = "predicate", indices = integer(0), ids = character(0),
        chromId = character(0), start = numeric(0), end = numeric(0),
        quo = quo, negate = FALSE
    )
}


# /// Negation Methods ///////////////////////////////////////////////

## ----
#' @rdname TaxaSelector-class
#' @aliases !,TaxaSelector-method
setMethod("!", "TaxaSelector", function(x) {
    x@negate <- !x@negate
    x
})

## ----
#' @rdname SiteSelector-class
#' @aliases !,SiteSelector-method
setMethod("!", "SiteSelector", function(x) {
    x@negate <- !x@negate
    x
})
