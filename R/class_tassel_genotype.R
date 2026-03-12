## ----
#' @title TasselGenotype Class
#' @description An S4 class to represent a Tassel Genotype object.
#'
#' @slot jRefObj
#' A reference to a Java object (`jobjRef`) associated with the genotype.
#' @slot jMemAddress
#' A character string representing the memory address of the Java object.
#' @slot jClass
#' A character string representing the Java class of the object.
#'
#' @details
#' This class is designed to interface with TASSEL 5 for genotype
#' data management and analysis. It provides a structure to store and
#' interact with Java objects used in TASSEL 5.
#'
#' @name TasselGenotype-class
#' @rdname TasselGenotype-class
#' @exportClass TasselGenotype
setClass(
    Class = "TasselGenotype",
    slots = c(
        jRefObj     = "jobjRef",
        jMemAddress = "character",
        jClass      = "character"
    )
)


## ----
#' @title
#' Read Genotype Data
#'
#' @description
#' This function reads genotype data from a file path or an R matrix.
#' It supports optional sorting of positions and retaining depth
#' information.
#'
#' @details
#' \itemize{
#'   \item If \code{x} is a character string:
#'     \itemize{
#'       \item The function checks if the file exists.
#'       \item Reads the genotype data from the file path using
#'       \code{readGenotypeFromPath}.
#'     }
#'   \item If \code{x} is a matrix:
#'     \itemize{
#'       \item The function processes the genotype data using
#'       \code{readGenotypeFromRMatrix}.
#'     }
#'   \item If \code{x} is neither a character string nor a matrix:
#'     \itemize{
#'       \item An error is raised.
#'     }
#' }
#'
#' @param x
#' A character string representing the file path to the genotype data
#' or a matrix containing genotype data.
#' @param sortPositions
#' A logical value indicating whether to sort positions in the
#' genotype data. Default is \code{FALSE}.
#' @param keepDepth
#' A logical value indicating whether to retain depth information in
#' the genotype data. Default is \code{FALSE}.
#'
#' @examples
#' \dontrun{
#' # Read genotype data from a file
#' readGenotype("path/to/genotype/file.txt", sortPositions = TRUE, keepDepth = TRUE)
#'
#' # Read genotype data from a matrix
#' genotypeMatrix <- matrix(data = ..., nrow = ..., ncol = ...)
#' readGenotype(genotypeMatrix)
#' }
#'
#' @return
#' A processed genotype object based on the input data.
#'
#' @export
readGenotype <- function(x, sortPositions = FALSE, keepDepth = FALSE) {
    if (is.character(x)) {
        xNorm <- normalizePath(x, mustWork = FALSE)
        if (!file.exists(xNorm)) {
            rlang::abort("The input path is not a valid file")
        }

        readGenotypeFromPath(xNorm, sortPositions, keepDepth)
    } else if (is.matrix(x)) {
        readNumericGenotypeFromRMatrix(x, asTGP = FALSE)
    } else {
        rlang::abort("Unsupported data type")
    }
}



# /// Methods (show) ////////////////////////////////////////////////

## ----
#' @title
#' Display TasselGenotype Object
#'
#' @description
#' This method is used to display a summary of a `TasselGenotype`
#' object. It prints genotype display information, including the
#' number of taxa, number of sites, and memory address of the Java
#' object.
#'
#' @param object
#' An object of class `TasselGenotype`.
#'
#' @method show TasselGenotype
setMethod("show", "TasselGenotype", function(object) {
    fgs <- formatGtStrings(object@jRefObj)
    printGtDisp(
        fgs    = fgs,
        nTaxa  = object@jRefObj$numberOfTaxa(),
        nSites = object@jRefObj$numberOfSites(),
        jMem   = object@jMemAddress
    )
})



# /// Methods (general) /////////////////////////////////////////////

## ----
#' @rdname javaRefObj
#' @export
setMethod(
    f = "javaRefObj",
    signature = signature(object = "TasselGenotype"),
    definition = function(object) {
        return(object@jRefObj)
    }
)



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


# /// Internal Helpers (bracket filtering) ///////////////////////////

## ----
# Wrap a filtered Java GenotypeTable back into the same R class as
# the original object.
newTasselGenotype <- function(jGt, original) {
    jMemAddress <- gsub(".*@", "", rJava::.jstrVal(jGt))
    jClass <- rJava::.jclass(jGt)
    methods::new(
        Class       = is(original)[[1]],
        jRefObj     = jGt,
        jMemAddress = jMemAddress,
        jClass      = jClass
    )
}


## ----
# Build a site metadata data frame from a Java GenotypeTable for
# predicate evaluation in sitesWhere().
buildSiteMetadata <- function(jGt) {
    nSites    <- jGt$numberOfSites()
    positions <- jGt$positions()
    idx       <- seq_len(nSites) - 1L

    jRc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    plArrays <- jRc$genotypeTableToPositionListOfArrays(positions)

    siteId <- vapply(idx, function(i) {
        jGt$siteName(as.integer(i))
    }, FUN.VALUE = character(1))

    mafVals <- vapply(idx, function(i) {
        jGt$minorAlleleFrequency(as.integer(i))
    }, FUN.VALUE = numeric(1))

    alleleCounts <- vapply(idx, function(i) {
        jGt$totalGametesNonMissingForSite(as.integer(i))
    }, FUN.VALUE = integer(1))

    nTaxa <- jGt$numberOfTaxa()
    hetVals <- vapply(idx, function(i) {
        jGt$heterozygousCount(as.integer(i)) / nTaxa
    }, FUN.VALUE = numeric(1))

    isIndelVals <- vapply(idx, function(i) {
        jGt$isIndel(as.integer(i))
    }, FUN.VALUE = logical(1))

    isBiallelicVals <- vapply(idx, function(i) {
        alleles <- jGt$alleles(as.integer(i))
        length(alleles) <= 2L
    }, FUN.VALUE = logical(1))

    data.frame(
        siteIndex    = idx,
        siteId       = siteId,
        chrom        = plArrays$chromosomes,
        pos          = as.numeric(plArrays$startPos),
        maf          = mafVals,
        alleleCount  = alleleCounts,
        het          = hetVals,
        isIndel      = isIndelVals,
        isBiallelic  = isBiallelicVals,
        stringsAsFactors = FALSE
    )
}


## ----
# Resolve a taxa selector (character vector or TaxaSelector) to a
# character vector of taxa IDs.
resolveTaxaIds <- function(jGt, selector) {
    if (is.character(selector)) return(selector)

    if (!methods::is(selector, "TaxaSelector")) {
        rlang::abort("Taxa selector must be a character vector or TaxaSelector")
    }

    if (selector@type == "ids") return(selector@ids)

    if (selector@type == "predicate") {
        nTaxa <- jGt$numberOfTaxa()
        taxaIds <- vapply(seq_len(nTaxa) - 1L, function(i) {
            jGt$taxaName(as.integer(i))
        }, FUN.VALUE = character(1))

        mask <- rlang::eval_tidy(
            selector@quo,
            data = data.frame(taxaId = taxaIds, stringsAsFactors = FALSE)
        )
        if (!is.logical(mask)) {
            rlang::abort("taxaWhere() expression must evaluate to a logical vector")
        }
        return(taxaIds[which(mask)])
    }

    rlang::abort(paste0("Unknown TaxaSelector type: ", selector@type))
}


## ----
# Apply a taxa selector to a Java GenotypeTable. Returns a filtered
# Java GenotypeTable.
applyTaxaSelector <- function(jGt, selector) {
    taxaIds <- resolveTaxaIds(jGt, selector)

    allTaxa <- vapply(seq_len(jGt$numberOfTaxa()) - 1L, function(i) {
        jGt$taxaName(as.integer(i))
    }, FUN.VALUE = character(1))

    negate <- methods::is(selector, "TaxaSelector") && selector@negate
    if (negate) {
        taxaIds <- setdiff(allTaxa, taxaIds)
    } else {
        missing <- setdiff(taxaIds, allTaxa)
        if (length(missing) == length(taxaIds)) {
            rlang::abort("No taxa match the selection criteria")
        }
        taxaIds <- intersect(taxaIds, allTaxa)
    }

    if (length(taxaIds) == 0) {
        rlang::abort("No taxa match the selection criteria")
    }

    builder <- rJava::.jnew("net.maizegenetics.taxa.TaxaListBuilder")
    builder$addAll(rJava::.jarray(taxaIds))
    taxaList <- builder$build()

    rJava::J("net.maizegenetics.dna.snp.FilterGenotypeTable")$getInstance(
        jGt, taxaList
    )
}


## ----
# Resolve a site selector to 0-based integer site indices.
resolveSiteIndices <- function(jGt, selector) {
    if (is.numeric(selector) || is.integer(selector)) {
        return(as.integer(selector))
    }

    if (is.character(selector)) {
        nSites <- jGt$numberOfSites()
        allNames <- vapply(seq_len(nSites) - 1L, function(i) {
            jGt$siteName(as.integer(i))
        }, FUN.VALUE = character(1))
        idx <- match(selector, allNames)
        return(as.integer(idx[!is.na(idx)] - 1L))
    }

    if (!methods::is(selector, "SiteSelector")) {
        rlang::abort("Site selector must be numeric, character, or SiteSelector")
    }

    switch(selector@type,
        "indices" = {
            selector@indices
        },
        "names" = {
            nSites <- jGt$numberOfSites()
            allNames <- vapply(seq_len(nSites) - 1L, function(i) {
                jGt$siteName(as.integer(i))
            }, FUN.VALUE = character(1))
            idx <- match(selector@ids, allNames)
            as.integer(idx[!is.na(idx)] - 1L)
        },
        "chrom" = {
            indices <- integer(0)
            chroms <- rJava::.jevalArray(jGt$chromosomes())
            for (jChrom in chroms) {
                if (jChrom$getName() %in% selector@chromId) {
                    rng <- jGt$firstLastSiteOfChromosome(jChrom)
                    indices <- c(indices, seq.int(rng[1], rng[2]))
                }
            }
            as.integer(indices)
        },
        "region" = {
            chroms <- rJava::.jevalArray(jGt$chromosomes())
            targetChrom <- NULL
            for (jChrom in chroms) {
                if (jChrom$getName() == selector@chromId) {
                    targetChrom <- jChrom
                    break
                }
            }
            if (is.null(targetChrom)) {
                rlang::abort(paste0("Chromosome '", selector@chromId, "' not found"))
            }

            chromRange <- jGt$firstLastSiteOfChromosome(targetChrom)

            posStart <- jGt$siteOfPhysicalPosition(
                as.integer(selector@start), targetChrom
            )
            if (posStart < 0) posStart <- -(posStart + 1L)

            posEnd <- jGt$siteOfPhysicalPosition(
                as.integer(selector@end), targetChrom
            )
            if (posEnd < 0) posEnd <- -(posEnd + 1L) - 1L

            posStart <- max(posStart, chromRange[1])
            posEnd   <- min(posEnd, chromRange[2])

            if (posStart > posEnd) return(integer(0))
            as.integer(seq.int(posStart, posEnd))
        },
        "predicate" = {
            meta <- buildSiteMetadata(jGt)
            mask <- rlang::eval_tidy(selector@quo, data = meta)
            if (!is.logical(mask)) {
                rlang::abort("sitesWhere() expression must evaluate to a logical vector")
            }
            as.integer(meta$siteIndex[which(mask)])
        },
        rlang::abort(paste0("Unknown SiteSelector type: ", selector@type))
    )
}


## ----
# Apply a site selector to a Java GenotypeTable. Returns a filtered
# Java GenotypeTable.
applySiteSelector <- function(jGt, selector) {
    indices <- resolveSiteIndices(jGt, selector)

    negate <- methods::is(selector, "SiteSelector") && selector@negate
    if (negate) {
        allIndices <- seq_len(jGt$numberOfSites()) - 1L
        indices <- as.integer(setdiff(allIndices, indices))
    }

    if (length(indices) == 0) {
        rlang::abort("No sites match the selection criteria")
    }

    rJava::J("net.maizegenetics.dna.snp.FilterGenotypeTable")$getInstance(
        jGt, rJava::.jarray(as.integer(indices))
    )
}


# /// Bracket Method /////////////////////////////////////////////////

## ----
#' @title Subset a TasselGenotype
#'
#' @description
#' Matrix-style subsetting for \code{TasselGenotype} objects using
#' the \code{gt[taxa, sites]} syntax.
#'
#' @param x A \code{TasselGenotype} object.
#' @param i Taxa selector: a character vector of IDs, a
#'   \code{\linkS4class{TaxaSelector}}, or missing.
#' @param j Site selector: an integer vector of 0-based indices, a
#'   character vector of site names, a
#'   \code{\linkS4class{SiteSelector}}, or missing.
#' @param ... Ignored.
#' @param drop Ignored.
#'
#' @return A new \code{TasselGenotype} (or subclass) containing the
#'   selected taxa and/or sites.
#'
#' @examples
#' \dontrun{
#' gt[taxa("B73", "Mo17"), ]
#' gt[, sites(0:999)]
#' gt[, sitesWhere(maf >= 0.05)]
#' gt[taxa("B73"), region("chr1", 1e6, 2e6)]
#' gt[, !sites(0:9)]
#' }
#'
#' @rdname TasselGenotype-class
#' @aliases [,TasselGenotype,ANY,ANY-method
setMethod("[", "TasselGenotype", function(x, i, j, ..., drop = FALSE) {
    jGt <- x@jRefObj
    if (!missing(i)) jGt <- applyTaxaSelector(jGt, i)
    if (!missing(j)) jGt <- applySiteSelector(jGt, j)
    newTasselGenotype(jGt, x)
})


