# /// S4 - class definition /////////////////////////////////////////

## ----
#' @title
#' TasselGenomicDataset Class Definition
#'
#' @description
#' Defines the \code{TasselGenomicDataset} class, which pairs genotype and
#' phenotype data that have been joined on their shared taxa. This is the
#' primary object for analyses that need both data types (e.g.
#' \code{assocModelFitter()}) and wraps TASSEL 5's \code{GenotypePhenotype}
#' class.
#'
#' @details
#' The \code{genotype} and \code{phenotype} slots hold wrappers around the
#' \emph{joined} tables rather than the tables originally handed to
#' \code{readGenomicDataset()}, so component accessors always reflect what the
#' dataset actually contains.
#'
#' @slot genotype
#' A \code{\linkS4class{TasselGenotype}} (or
#' \code{\linkS4class{TasselNumericGenotype}}) holding the joined genotype
#' table.
#' @slot phenotype
#' A \code{\linkS4class{TasselPhenotype}} holding the joined phenotype table.
#' @slot jRefObj
#' A \code{jobjRef} referencing the TASSEL 5 \code{GenotypePhenotype} object.
#' @slot jMemAddress
#' A \code{character} string representing the memory address of the Java
#' object.
#' @slot jClass
#' A \code{character} string representing the Java class name of the object.
#'
#' @include class_tassel_genotype.R
#' @include class_tassel_phenotype.R
#'
#' @name TasselGenomicDataset-class
#' @rdname TasselGenomicDataset-class
#' @exportClass TasselGenomicDataset
setClass(
    Class = "TasselGenomicDataset",
    slots = c(
        genotype    = "TasselGenotype",
        phenotype   = "TasselPhenotype",
        jRefObj     = "jobjRef",
        jMemAddress = "character",
        jClass      = "character"
    )
)



# /// Constructors //////////////////////////////////////////////////

## ----
# Create a TasselGenomicDataset from a Java GenotypePhenotype object
#
# @description
# This internal function builds a `TasselGenomicDataset` S4 object from a
# Java `GenotypePhenotype` object, wrapping the joined genotype and
# phenotype tables it carries.
#
# @param jGp
# A Java `GenotypePhenotype` object reference.
#
# @return
# An S4 object of class `TasselGenomicDataset`.
createTasselGenomicDataset <- function(jGp) {
    if (!methods::is(jGp, "jobjRef") || rJava::is.jnull(jGp)) {
        rlang::abort("`jGp` must be a non-null Java object reference")
    }

    jGt <- getGenotypeTable(jGp)
    jPh <- getPhenotypeTable(jGp)

    if (rJava::is.jnull(jGt) || rJava::is.jnull(jPh)) {
        rlang::abort(
            "`jGp` must hold both genotype and phenotype data"
        )
    }

    methods::new(
        Class       = "TasselGenomicDataset",
        genotype    = createTasselGenotype(jGt),
        phenotype   = createTasselPhenotype(jPh),
        jRefObj     = jGp,
        jMemAddress = gsub(".*@", "", rJava::.jstrVal(jGp)),
        jClass      = rJava::.jclass(jGp)
    )
}


## ----
# Join a Java GenotypeTable and Phenotype into a GenotypePhenotype
#
# @param jGt
# A Java `GenotypeTable` object reference.
# @param jPh
# A Java `Phenotype` object reference.
# @param join
# Either `"intersect"` (taxa present in both tables) or `"union"` (taxa
# present in either table).
#
# @return
# A Java `GenotypePhenotype` object reference.
joinGenotypePhenotype <- function(jGt, jPh, join = "intersect") {
    builder <- rJava::.jnew(TASSEL_JVM$GENO_PHENO_BUILDER)$
        genotype(jGt)$
        phenotype(jPh)

    builder <- if (join == "union") builder$union() else builder$intersect()

    # A Java 'Throwable' cannot be used as 'parent': rJava's '$' method
    # intercepts the fields rlang probes when formatting a chained error
    joinFailed <- function(cnd = NULL) {
        rlang::abort(c(
            "Could not join the genotype and phenotype data",
            "i" = "Do the two data sets share any taxa IDs?",
            if (!is.null(cnd)) {
                c("i" = paste0("TASSEL reported: ", conditionMessage(cnd)))
            }
        ))
    }

    # A join that leaves no taxa behind surfaces as a Java exception rather
    # than an empty result, so both outcomes report the likely cause
    jGp <- tryCatch(builder$build(), error = joinFailed)

    if (rJava::is.jnull(jGp)) {
        joinFailed()
    }

    return(jGp)
}


## ----
# Coerce the 'genotype' argument of readGenomicDataset() to a Java table
#
# @param x
# A `TasselGenotype` object, a path to a genotype file, or an R matrix of
# numeric genotype values.
#
# @return
# A Java `GenotypeTable` object reference.
resolveGenotypeArg <- function(x) {
    if (is.character(x) || is.matrix(x)) {
        x <- readGenotype(x)
    }

    jGt <- getGenotypeTable(x)

    if (rJava::is.jnull(jGt)) {
        rlang::abort(c(
            "`genotype` does not contain genotype data",
            "i" = paste0(
                "Pass a <TasselGenotype> object, a path to a genotype file, ",
                "or a numeric matrix."
            )
        ))
    }

    return(jGt)
}


## ----
# Coerce the 'phenotype' argument of readGenomicDataset() to a Java table
#
# @param x
# A `TasselPhenotype` object, a path to a phenotype file, or a `data.frame`.
# @param attr
# Attribute metadata, required when `x` is a `data.frame`. See
# `readPhenotype()`.
#
# @return
# A Java `Phenotype` object reference.
resolvePhenotypeArg <- function(x, attr = NULL) {
    if (is.character(x) || is.data.frame(x)) {
        x <- readPhenotype(x, attr = attr)
    }

    jPh <- getPhenotypeTable(x)

    if (rJava::is.jnull(jPh)) {
        rlang::abort(c(
            "`phenotype` does not contain phenotype data",
            "i" = paste0(
                "Pass a <TasselPhenotype> object, a path to a phenotype ",
                "file, or a `data.frame` along with `attr`."
            )
        ))
    }

    return(jPh)
}


## ----
#' @title
#' Read genotype and phenotype data into a single joined dataset
#'
#' @description
#' Joins genotype and phenotype data on their shared taxa and returns a
#' \code{\linkS4class{TasselGenomicDataset}} object. Both arguments accept
#' either already-constructed \code{rTASSEL} objects or the raw inputs that
#' \code{\link{readGenotype}()} and \code{\link{readPhenotype}()} understand,
#' so a dataset can be built in one call.
#'
#' @param genotype
#' A \code{\linkS4class{TasselGenotype}} object, a path to a genotype file
#' (e.g. \code{*.vcf}, \code{*.hmp.txt}), or a numeric matrix of genotype
#' values. Paths and matrices are passed to \code{\link{readGenotype}()}.
#' @param phenotype
#' A \code{\linkS4class{TasselPhenotype}} object, a path to a phenotype file,
#' or a \code{data.frame} of phenotype data. Paths and data frames are passed
#' to \code{\link{readPhenotype}()}.
#' @param attr
#' Attribute metadata describing the columns of \code{phenotype}. Required
#' when \code{phenotype} is a \code{data.frame} and ignored otherwise. See
#' \code{\link{readPhenotype}()} for the expected format.
#' @param join
#' How taxa from the two data sets are combined. \code{"intersect"} (the
#' default) keeps only taxa found in both; \code{"union"} keeps taxa found in
#' either.
#'
#' @return
#' A \code{\linkS4class{TasselGenomicDataset}} object.
#'
#' @examples
#' \dontrun{
#' # From file paths
#' ds <- readGenomicDataset(
#'     genotype  = "path/to/genotype.hmp.txt",
#'     phenotype = "path/to/phenotype.txt"
#' )
#'
#' # From existing rTASSEL objects
#' gt <- readGenotype("path/to/genotype.hmp.txt")
#' ph <- readPhenotype("path/to/phenotype.txt")
#' ds <- readGenomicDataset(gt, ph)
#'
#' # Keep taxa found in either data set
#' ds <- readGenomicDataset(gt, ph, join = "union")
#' }
#'
#' @export
readGenomicDataset <- function(
    genotype,
    phenotype,
    attr = NULL,
    join = c("intersect", "union")
) {
    join <- rlang::arg_match(join)

    jGt <- resolveGenotypeArg(genotype)
    jPh <- resolvePhenotypeArg(phenotype, attr)

    createTasselGenomicDataset(joinGenotypePhenotype(jGt, jPh, join))
}



# /// Methods (show) ////////////////////////////////////////////////

## ----
# Collapse a phenotype attribute summary into a single display string
#
# @param attrSummary
# The `attrSummary` slot of a `TasselPhenotype`: a named list of trait
# counts keyed by TASSEL attribute type.
#
# @return
# A single `character` value, e.g. `"data: 3, taxa: 1"`.
formatTraitSummary <- function(attrSummary) {
    if (length(attrSummary) == 0) {
        return("no traits")
    }

    paste0(names(attrSummary), ": ", unlist(attrSummary), collapse = ", ")
}


## ----
#' @title
#' Display summary information of a TasselGenomicDataset object
#'
#' @description
#' Prints the dimensions of the joined genotype table, a breakdown of the
#' phenotype traits by TASSEL attribute type, and the memory address of the
#' backing Java object.
#'
#' @param object
#' A \code{TasselGenomicDataset} object.
#'
#' @rdname TasselGenomicDataset-class
#' @aliases show,TasselGenomicDataset-method
setMethod("show", "TasselGenomicDataset", function(object) {
    jGt <- object@genotype@jRefObj

    infoLine <- function(...) {
        pillar::style_subtle(paste0("# ", cli::symbol$info, " ", ...))
    }

    cli::cat_line(pillar::style_subtle(paste0(
        "# A ", cli::style_bold("TasselGenomicDataset"), " object: ",
        jGt$numberOfTaxa(), " taxa ", cli::symbol$times, " ",
        jGt$numberOfSites(), " sites"
    )))
    cli::cat_line()
    cli::cat_line(infoLine(
        "Genotype..: <", methods::is(object@genotype)[[1]], ">"
    ))
    cli::cat_line(infoLine(
        "Phenotype.: ", nrow(object@phenotype@attrData), " traits (",
        formatTraitSummary(object@phenotype@attrSummary), ")"
    ))
    cli::cat_line(infoLine(
        "Java memory address: 0x", cli::style_bold(object@jMemAddress)
    ))
})



# /// Methods (component accessors) /////////////////////////////////

## ----
#' @rdname genotype
#' @aliases genotype,TasselGenomicDataset-method
#' @export
setMethod(
    f = "genotype",
    signature = signature(object = "TasselGenomicDataset"),
    definition = function(object) {
        return(object@genotype)
    }
)


## ----
#' @rdname phenotype
#' @aliases phenotype,TasselGenomicDataset-method
#' @export
setMethod(
    f = "phenotype",
    signature = signature(object = "TasselGenomicDataset"),
    definition = function(object) {
        return(object@phenotype)
    }
)


## ----
#' @rdname javaRefObj
#' @aliases javaRefObj,TasselGenomicDataset-method
#' @export
setMethod(
    f = "javaRefObj",
    signature = signature(object = "TasselGenomicDataset"),
    definition = function(object) {
        return(object@jRefObj)
    }
)



# /// Methods (taxa + positions) ////////////////////////////////////

## ----
#' @rdname taxaList
#' @aliases taxaList,TasselGenomicDataset-method
#' @export
setMethod("taxaList", "TasselGenomicDataset", function(tasObj) {
    taxaList(genotype(tasObj))
})


## ----
#' @rdname positionList
#' @aliases positionList,TasselGenomicDataset-method
#' @export
setMethod("positionList", "TasselGenomicDataset", function(tasObj) {
    positionList(genotype(tasObj))
})


## ----
#' @rdname seqnames
#' @aliases seqnames,TasselGenomicDataset-method
#' @export
setMethod("seqnames", "TasselGenomicDataset", function(x) {
    seqnames(genotype(x))
})



# /// Methods (summary) /////////////////////////////////////////////

## ----
#' @rdname siteSummary
#' @aliases siteSummary,TasselGenomicDataset-method
#' @export
setMethod("siteSummary", "TasselGenomicDataset", function(tasObj) {
    siteSummary(genotype(tasObj))
})


## ----
#' @rdname taxaSummary
#' @aliases taxaSummary,TasselGenomicDataset-method
#' @export
setMethod("taxaSummary", "TasselGenomicDataset", function(tasObj) {
    taxaSummary(genotype(tasObj))
})


## ----
#' @rdname attributeData
#' @aliases attributeData,TasselGenomicDataset-method
#' @export
setMethod(
    f = "attributeData",
    signature = signature(object = "TasselGenomicDataset"),
    definition = function(object) {
        attributeData(phenotype(object))
    }
)


## ----
#' @rdname traitNames
#' @aliases traitNames,TasselGenomicDataset-method
#' @export
setMethod(
    f = "traitNames",
    signature = signature(object = "TasselGenomicDataset"),
    definition = function(object) {
        traitNames(phenotype(object))
    }
)



# /// Methods (coercion) ////////////////////////////////////////////

## ----
#' @title Coerce a genomic dataset's phenotype data to a data frame
#'
#' @description
#' Returns the joined phenotype data held by a
#' \code{\linkS4class{TasselGenomicDataset}} as a \code{tibble}.
#'
#' @param x A \code{TasselGenomicDataset} object.
#' @param row.names Ignored, present for generic compatibility.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Additional arguments to be passed to or from methods.
#'
#' @return A \code{tibble} of phenotype data.
#'
#' @export
as.data.frame.TasselGenomicDataset <- function(
    x,
    row.names = NULL,
    optional = FALSE,
    ...
) {
    as.data.frame(phenotype(x), row.names = row.names, optional = optional, ...)
}


## ----
#' @title Coerce a genomic dataset's genotype data to an R matrix
#'
#' @description
#' Returns the joined genotype table held by a
#' \code{\linkS4class{TasselGenomicDataset}} as a matrix of dosage values.
#'
#' @param x A \code{TasselGenomicDataset} object.
#' @param ... Additional arguments to be passed to or from methods.
#'
#' @return An \code{integer} matrix of taxa (rows) by sites (columns).
#'
#' @export
as.matrix.TasselGenomicDataset <- function(x, ...) {
    as.matrix(genotype(x), ...)
}



# /// Bracket Method /////////////////////////////////////////////////

## ----
#' @title Subset a TasselGenomicDataset
#'
#' @description
#' Matrix-style subsetting for \code{TasselGenomicDataset} objects using the
#' \code{ds[taxa, sites]} syntax. Selectors are applied to the genotype table
#' and the phenotype data is then re-joined against the surviving taxa, so
#' both components of the returned dataset stay in step.
#'
#' @param x A \code{TasselGenomicDataset} object.
#' @param i Taxa selector: a character vector of IDs, a
#'   \code{\linkS4class{TaxaSelector}}, or missing.
#' @param j Site selector: an integer vector of 1-based indices, a
#'   character vector of site names, a
#'   \code{\linkS4class{SiteSelector}}, or missing.
#' @param ... Ignored.
#' @param drop Ignored.
#'
#' @return A new \code{TasselGenomicDataset} containing the selected taxa
#'   and/or sites.
#'
#' @examples
#' \dontrun{
#' ds[taxa("B73", "Mo17"), ]
#' ds[, sitesWhere(maf >= 0.05)]
#' ds[taxaWhere(startsWith(taxaId, "NAM")), region("chr1", 1e6, 2e6)]
#' }
#'
#' @rdname TasselGenomicDataset-class
#' @aliases [,TasselGenomicDataset,ANY,ANY-method
setMethod("[", "TasselGenomicDataset", function(x, i, j, ..., drop = FALSE) {
    jGt <- x@genotype@jRefObj
    if (!missing(i)) jGt <- applyTaxaSelector(jGt, i)
    if (!missing(j)) jGt <- applySiteSelector(jGt, j)

    createTasselGenomicDataset(
        joinGenotypePhenotype(jGt, x@phenotype@jRefObj)
    )
})
