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



# /// Methods (taxa + positions) ////////////////////////////////////

## ----
#' @rdname taxaList
#' @aliases taxaList,TasselGenotype-method
#' @export
setMethod("taxaList", "TasselGenotype", function(tasObj) {
    .taxaNames(tasObj@jRefObj$taxa())
})


## ----
#' @rdname positionList
#' @aliases positionList,TasselGenotype-method
#' @export
setMethod("positionList", "TasselGenotype", function(tasObj) {
    sites <- rJava::new(
        rJava::J("net.maizegenetics.dna.map.PositionListTableReport"),
        tasObj@jRefObj$positions()
    )
    tableReportToDF(sites)
})


## ----
#' @title Chromosome (sequence) IDs from TasselGenotype
#'
#' @description
#' Returns the chromosome (sequence) IDs present in the genotype table,
#' in the order returned by TASSEL.
#'
#' @param x A \code{TasselGenotype} object.
#'
#' @return A character vector of chromosome IDs.
#'
#' @rdname seqnames
#' @aliases seqnames,TasselGenotype-method
#' @export
setMethod("seqnames", "TasselGenotype", function(x) {
    .jStrings(javaRefObj(x)$chromosomes())
})


## ----
#' @title Marker positions as a GRanges object
#'
#' @description
#' Returns the positions of a genotype table as a
#' \code{GenomicRanges::GRanges} object, which is the form the rest of
#' Bioconductor expects. Each marker is a width-1 range.
#'
#' @details
#' The result can be handed straight back to \code{\link{region}()} or
#' \code{\link{overlaps}()} to filter on, and carries the same
#' information as \code{\link{positionList}()} in a different shape.
#'
#' @param x A \code{TasselGenotype}, \code{TasselGenomicDataset}, or
#'   deprecated \code{TasselGenotypePhenotype} object.
#' @param use.names Name each range with its marker ID? Defaults to
#'   \code{TRUE}.
#' @param use.mcols Carry the remaining \code{\link{positionList}()}
#'   columns (\code{Site}, \code{Name}, and \code{VARIANT}) across as
#'   metadata columns? Defaults to \code{FALSE}.
#' @param ... Additional arguments, for use in specific methods.
#'
#' @return A \code{GRanges} object with one range per marker.
#'
#' @examples
#' \dontrun{
#' granges(gt)
#'
#' granges(gt, use.mcols = TRUE)
#'
#' # Filter one genotype table on the ranges of another
#' gt[, sitesWhere(overlaps(granges(otherGt)))]
#' }
#'
#' @rdname granges
#' @aliases granges,TasselGenotype-method
#' @export
setMethod(
    "granges",
    "TasselGenotype",
    function(x, use.names = TRUE, use.mcols = FALSE, ...) {
        .positionRanges(positionList(x), use.names, use.mcols)
    }
)



# /// Methods (summary) /////////////////////////////////////////////

## ----
#' @rdname siteSummary
#' @aliases siteSummary,TasselGenotype-method
#' @export
setMethod("siteSummary", "TasselGenotype", function(tasObj) {
    .runGenotypeSummary(tasObj@jRefObj, doSite = TRUE)
})


## ----
#' @rdname taxaSummary
#' @aliases taxaSummary,TasselGenotype-method
#' @export
setMethod("taxaSummary", "TasselGenotype", function(tasObj) {
    .runGenotypeSummary(tasObj@jRefObj, doTaxa = TRUE)
})



# /// Methods (coercion) ////////////////////////////////////////////

## ----
#' @title Coerce genotype data to an R matrix
#'
#' @description
#' Converts the genotype table held by a \code{TasselGenotype} object into a
#' matrix with taxa as rows and sites as columns.
#'
#' @details
#' Two readings of the same calls are available. \code{"dosage"} counts the
#' alternate alleles a taxon carries at a site, which is the form most
#' models want. \code{"allele"} returns the calls as TASSEL itself spells
#' them, which for nucleotide data means one IUPAC code per cell and
#' \code{"N"} for a missing call.
#'
#' Both readings are as large as the data, so a table of any real size
#' should be filtered down to the taxa and sites of interest first.
#'
#' @param x A \code{TasselGenotype} object.
#' @param type Reading of the genotype calls to return. Either
#'   \code{"dosage"} (the default) for alternate allele counts, or
#'   \code{"allele"} for call strings.
#' @param ... Additional arguments to be passed to or from methods.
#'
#' @return
#' An \code{integer} matrix of taxa (rows) by sites (columns) when
#' \code{type} is \code{"dosage"}, or a \code{character} matrix of the same
#' shape when \code{type} is \code{"allele"}.
#'
#' @examples
#' \dontrun{
#' as.matrix(gt)
#' as.matrix(gt, type = "allele")
#' }
#'
#' @export
as.matrix.TasselGenotype <- function(x, type = c("dosage", "allele"), ...) {
    type <- match.arg(type)

    if (!x@jRefObj$hasGenotype()) {
        rlang::abort(c(
            "`x` does not contain discrete genotype calls",
            "i" = "Only allele-based genotype tables can be coerced to a matrix"
        ))
    }

    taxa <- taxaList(x)
    siteNames <- positionList(x)$Name

    switch(
        type,
        "dosage" = .dosageMatrix(x@jRefObj, taxa = taxa, siteNames = siteNames),
        "allele" = .alleleStringMatrix(x@jRefObj, taxa = taxa, siteNames = siteNames)
    )
}


## ----
#' @title Coerce genotype data to a SummarizedExperiment
#'
#' @description
#' Assembles the dosage matrix, marker positions, and taxa IDs of a
#' \code{TasselGenotype} into a
#' \code{SummarizedExperiment::SummarizedExperiment}, which is the
#' container the rest of Bioconductor expects.
#'
#' @details
#' A \code{SummarizedExperiment} puts features in rows and samples in
#' columns, so the assay is the transpose of
#' \code{\link{as.matrix.TasselGenotype}}: sites are rows and taxa are
#' columns.
#'
#' This replaces the deprecated \code{\link{getSumExpFromGenotypeTable}()}.
#'
#' @param from A \code{TasselGenotype} object.
#' @param to The target class, \code{"SummarizedExperiment"}.
#' @param strict Supplied by \code{\link[methods]{as}()}; unused here.
#'
#' @return
#' A \code{SummarizedExperiment} of sites (rows) by taxa (columns).
#'
#' @name coerce-TasselGenotype-SummarizedExperiment
#' @aliases coerce,TasselGenotype,SummarizedExperiment-method
#'
#' @examples
#' \dontrun{
#' se <- as(gt, "SummarizedExperiment")
#' }
#'
#' @export
setAs("TasselGenotype", "SummarizedExperiment", function(from) {
    .genotypeSummarizedExperiment(from)
})



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
#' @param j Site selector: an integer vector of 1-based indices, a
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
#' gt[, sites(1:1000)]
#' gt[, sitesWhere(maf >= 0.05)]
#' gt[taxaWhere(notMissing >= 0.8), ]
#' gt[taxa("B73"), region("chr1", 1e6, 2e6)]
#' gt[, !sites(1:10)]
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
