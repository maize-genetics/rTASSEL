# /// JVM -> R transport helpers /////////////////////////////////////
#
# Every route data takes out of the JVM lives here, so that the rules
# each one follows - missing value codes, orientation, dimnames - are
# stated once and shared by every accessor and coercion method that
# needs them.


## ----
#' @title Character vector from a collection of Java objects
#'
#' @description
#' \code{rJava} hands back object arrays and collections as R lists of
#' \code{jobjRef}. This resolves such a list to the \code{toString()} of
#' each element.
#'
#' @param x A list or vector of \code{jobjRef} objects.
#'
#' @return A \code{character} vector the same length as \code{x}.
#'
#' @noRd
#' @importFrom rJava .jstrVal
.jStrings <- function(x) {
    vapply(x, rJava::.jstrVal, FUN.VALUE = character(1), USE.NAMES = FALSE)
}


## ----
#' @title Taxa IDs from a Java taxa list
#'
#' @description
#' The single route taxa IDs take into R, shared by the \code{taxaList()}
#' methods and the legacy \code{getTaxaIDs()} helper.
#'
#' @param jTaxaList A Java \code{TaxaList} reference.
#'
#' @return A \code{character} vector of taxa IDs.
#'
#' @noRd
#' @importFrom rJava J
.taxaNames <- function(jTaxaList) {
    rJava::J(TASSEL_JVM$R_METHODS)$genotypeTableToSampleNameArray(jTaxaList)
}


## ----
#' @title Dosage matrix from a Java genotype table
#'
#' @description
#' Returns TASSEL's dosage byte array as an R matrix of taxa (rows) by
#' sites (columns).
#'
#' @param jGt A Java \code{GenotypeTable} reference.
#' @param taxa A \code{character} vector of taxa IDs, or \code{NULL}.
#' @param siteNames A \code{character} vector of marker names, or
#'   \code{NULL}.
#' @param asInteger Coerce to \code{integer} and map TASSEL's missing
#'   value to \code{NA}? When \code{FALSE} the \code{raw} bytes are
#'   returned as they arrive, which uses a quarter of the memory but
#'   carries no \code{NA}.
#'
#' @return An \code{integer} (or \code{raw}) matrix.
#'
#' @noRd
#' @importFrom rJava .jevalArray
#' @importFrom rJava J
.dosageMatrix <- function(jGt, taxa = NULL, siteNames = NULL, asInteger = TRUE) {
    jrc <- rJava::J(TASSEL_JVM$R_METHODS)
    m <- rJava::.jevalArray(
        jrc$genotypeTableToDosageByteArray(jGt),
        simplify = TRUE
    )

    if (asInteger) {
        mode(m) <- "integer"

        # 128 is a conversion artifact of TASSEL's unsigned missing value
        m[m == 128] <- NA
    }

    if (!is.null(taxa) || !is.null(siteNames)) {
        dimnames(m) <- list(taxa, siteNames)
    }

    return(m)
}


## ----
#' @title Allele call strings from a Java genotype table
#'
#' @description
#' Returns the genotype calls as a \code{character} matrix of taxa (rows)
#' by sites (columns), in the same spelling TASSEL itself reports (IUPAC
#' codes for nucleotide data).
#'
#' @details
#' TASSEL packs a diploid call into a single byte, so a table holds at
#' most 256 distinct calls no matter how large it is. One call to
#' \code{genotypeAllSites()} per taxon brings the codes across, and each
#' distinct code is then spelled out once by asking the table to render a
#' cell that holds it.
#'
#' @param jGt A Java \code{GenotypeTable} reference.
#' @param taxa A \code{character} vector of taxa IDs, or \code{NULL}.
#' @param siteNames A \code{character} vector of marker names, or
#'   \code{NULL}.
#'
#' @return A \code{character} matrix.
#'
#' @noRd
.alleleStringMatrix <- function(jGt, taxa = NULL, siteNames = NULL) {
    nTaxa <- jGt$numberOfTaxa()
    nSites <- jGt$numberOfSites()

    codes <- matrix(
        unlist(
            lapply(
                seq_len(nTaxa),
                function(i) as.integer(jGt$genotypeAllSites(as.integer(i - 1L)))
            ),
            use.names = FALSE
        ),
        nrow = nSites,
        ncol = nTaxa
    )

    m <- .decodeGenotypeCodes(jGt, codes, nSites)
    dimnames(m) <- list(taxa, siteNames)

    return(m)
}


## ----
#' @title Spell out packed genotype byte codes
#'
#' @description
#' Maps a sites-by-taxa matrix of TASSEL's packed genotype bytes to the
#' call strings the table reports for them.
#'
#' @param jGt A Java \code{GenotypeTable} reference.
#' @param codes An \code{integer} matrix of sites (rows) by taxa
#'   (columns).
#' @param nSites Number of sites in \code{jGt}.
#'
#' @return A \code{character} matrix of taxa (rows) by sites (columns).
#'
#' @noRd
#' @importFrom rJava .jevalArray
.decodeGenotypeCodes <- function(jGt, codes, nSites) {
    # TASSEL reports a single row of allele states when every site shares
    # them, which is the case for all nucleotide data. Only when the
    # states vary from site to site does a code need decoding per site.
    sharedStates <- length(rJava::.jevalArray(jGt$alleleDefinitions())) == 1L

    decodeAt <- function(cellIdx) {
        taxonIdx <- (cellIdx - 1L) %/% nSites
        siteIdx <- (cellIdx - 1L) %% nSites
        jGt$genotypeAsString(as.integer(taxonIdx), as.integer(siteIdx))
    }

    if (sharedStates) {
        distinct <- sort(unique(as.vector(codes)))
        spelling <- vapply(
            match(distinct, codes),
            decodeAt,
            FUN.VALUE = character(1)
        )

        return(t(matrix(
            spelling[match(codes, distinct)],
            nrow = nSites
        )))
    }

    out <- matrix(NA_character_, nrow = nSites, ncol = ncol(codes))
    for (site in seq_len(nSites)) {
        siteCodes <- codes[site, ]
        distinct <- sort(unique(siteCodes))
        spelling <- vapply(
            distinct,
            function(cd) {
                jGt$genotypeAsString(
                    as.integer(match(cd, siteCodes) - 1L),
                    as.integer(site - 1L)
                )
            },
            FUN.VALUE = character(1)
        )
        out[site, ] <- spelling[match(siteCodes, distinct)]
    }

    return(t(out))
}


## ----
#' @title Reference probability matrix from a Java genotype table
#'
#' @description
#' Returns the reference probabilities of a numeric genotype table as an
#' R matrix of taxa (rows) by sites (columns).
#'
#' @details
#' TASSEL 5 exposes reference probabilities one cell at a time, so this
#' costs one JNI call per cell and is only practical on tables that have
#' been filtered down first.
#'
#' @param jGt A Java \code{GenotypeTable} reference.
#' @param taxa A \code{character} vector of taxa IDs, or \code{NULL}.
#' @param siteNames A \code{character} vector of marker names, or
#'   \code{NULL}.
#'
#' @return A \code{numeric} matrix.
#'
#' @noRd
.refProbMatrix <- function(jGt, taxa = NULL, siteNames = NULL) {
    nTaxa <- jGt$numberOfTaxa()
    nSites <- jGt$numberOfSites()

    m <- matrix(
        unlist(
            lapply(
                seq_len(nTaxa),
                function(i) {
                    vapply(
                        seq_len(nSites),
                        function(j) {
                            jGt$referenceProbability(
                                as.integer(i - 1L),
                                as.integer(j - 1L)
                            )
                        },
                        FUN.VALUE = numeric(1)
                    )
                }
            ),
            use.names = FALSE
        ),
        nrow = nTaxa,
        ncol = nSites,
        byrow = TRUE,
        dimnames = list(taxa, siteNames)
    )

    return(m)
}


## ----
#' @title Genomic ranges from a position table
#'
#' @description
#' Turns the table returned by \code{positionList()} into a
#' \code{GRanges} object of width-1 ranges, one per marker.
#'
#' @param pl A \code{tibble} as returned by \code{positionList()}.
#' @param use.names Name each range with its marker ID?
#' @param use.mcols Carry the remaining position columns across as
#'   metadata columns?
#'
#' @return A \code{GRanges} object.
#'
#' @noRd
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges granges
#' @importFrom IRanges IRanges
.positionRanges <- function(pl, use.names = TRUE, use.mcols = FALSE) {
    args <- list(
        seqnames = pl$Chromosome,
        ranges   = IRanges::IRanges(start = pl$Position, width = 1L)
    )

    if (use.mcols) {
        rangeCols <- c("Chromosome", "Position")
        args <- c(args, as.list(pl[setdiff(colnames(pl), rangeCols)]))
    }

    gr <- do.call(GenomicRanges::GRanges, args)

    if (use.names) {
        names(gr) <- pl$Name
    }

    return(gr)
}


## ----
#' @title SummarizedExperiment from genotype data
#'
#' @description
#' Assembles the dosage matrix, marker positions, and taxa IDs of a
#' genotype table into a \code{SummarizedExperiment}.
#'
#' @param x A \code{TasselGenotype} or \code{TasselGenomicDataset}
#'   object.
#'
#' @return A \code{SummarizedExperiment} object of sites (rows) by taxa
#'   (columns).
#'
#' @noRd
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
.genotypeSummarizedExperiment <- function(x) {
    taxa <- taxaList(x)

    SummarizedExperiment::SummarizedExperiment(
        # Features go in rows and samples in columns, which is the
        # transpose of the taxa-by-sites dosage matrix
        assays    = t(as.matrix(x)),
        rowRanges = granges(x, use.mcols = TRUE),
        colData   = S4Vectors::DataFrame(Sample = taxa, row.names = taxa)
    )
}


## ----
#' @title Pairwise matrix from a Java distance matrix
#'
#' @description
#' Returns a TASSEL \code{DistanceMatrix} as an ordinary \code{numeric}
#' matrix.
#'
#' @param jDist A Java \code{DistanceMatrix} reference.
#' @param taxa A \code{character} vector of taxa IDs, or \code{NULL}.
#'
#' @return A \code{numeric} matrix.
#'
#' @noRd
#' @importFrom rJava .jevalArray
.distanceToMatrix <- function(jDist, taxa = NULL) {
    m <- rJava::.jevalArray(jDist$getDistances(), simplify = TRUE)
    dimnames(m) <- list(taxa, taxa)

    return(m)
}
