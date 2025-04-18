# /// Display functions (core) //////////////////////////////////////

## ----
# Generate Ellipsis with Optional Spacing
#
# @description
# Creates an ellipsis symbol with optional spacing on either side.
#
# @param ind
# Integer indicating the number of spaces to add on each side. If
# \code{NULL}, no spacing is added.
#
# @return
# A character string containing a formatted ellipsis
subEllipsis <- function(ind = NULL) {
    spacer <- ifelse(!is.null(ind), strrep(" ", ind), "")
    sprintf(
        paste0(spacer, "%s", spacer),
        pillar::style_subtle(cli::symbol$ellipsis)
    )
}


## ----
# Truncate Genotype Name
#
# @description
# Truncates a genotype identifier if it exceeds the maximum length,
# adding an ellipsis.
#
# @param id
# Character string to be truncated
# @param maxLength
# Maximum allowed length for the string
#
# @return
# A character string, either truncated with ellipsis or unchanged
truncateGtName <- function(id, maxLength = 5) {
    if (nchar(id) > maxLength) {
        paste0(substr(id, 1, maxLength - 1), cli::symbol$ellipsis)
    } else {
        id
    }
}


## ----
# Format Genotype Names for Display
#
# @description
# Formats a list of genotype taxa names for display, handling
# truncation and padding. If there are more taxa than the display
# limit, shows first n taxa, ellipsis, and last taxon.
#
# @param gt
# Genotype object containing taxa names
# @param nTaxa
# Number of taxa to display before truncating
# @param maxLength
# Maximum length allowed for each taxa name
#
# @return
# A character vector of formatted and padded taxa names
formatGtNames <- function(gt, nTaxa = 5, maxLength = 5) {
    maxTaxa <- gt$numberOfTaxa()

    if (maxTaxa > nTaxa + 2) {
        taxaIds <- vapply(seq_len(nTaxa), function(it) {
            gt$taxaName(as.integer(it - 1))
        }, FUN.VALUE = character(1))

        taxaIds <- c(
            "",
            taxaIds,
            paste0(strrep(" ", maxLength - 1), cli::symbol$ellipsis),
            gt$taxaName(as.integer(gt$numberOfTaxa() - 1))
        )
    } else {
        taxaIds <- vapply(seq_len(maxTaxa), function(it) {
            gt$taxaName(as.integer(it - 1))
        }, FUN.VALUE = character(1))

        taxaIds <- c(
            "",
            taxaIds
        )
    }

    truncTaxaIds <- vapply(
        taxaIds,
        truncateGtName, maxLength = maxLength,
        FUN.VALUE = character(1)
    )

    paddedTrunTaxaIds <- pillar::style_subtle(
        sprintf(
            paste0("%", maxLength, "s "),
            truncTaxaIds
        )
    )

    return(paddedTrunTaxaIds)
}


## ----
# Add Row and Column IDs to Genotype Display
#
# @description
# Adds row and column identifiers to the genotype display matrix
#
# @param gt
# Genotype object
# @param fgs
# Formatted genotype strings
# @param nSites
# Number of sites to display before truncating
# @param numeric
# Logical indicating if numeric formatting should be used
#
# @return
# A vector of formatted strings with row and column IDs
addGtRowColIds <- function(gt, fgs, nSites = 10, numeric = FALSE) {
    maxSite <- gt$numberOfSites()
    fmt <- if (numeric) "%2d     " else "%2d "

    if (maxSite > nSites + 2) {
        colHeader <- c(
            sprintf(fmt, seq_len(nSites) - 1),
            sprintf(" %s ", cli::symbol$ellipsis),
            sprintf(" %d", maxSite - 1)
        )
    } else {
        colHeader <- sprintf(fmt, seq_len(maxSite) - 1)
    }

    colHeader <- pillar::style_subtle(paste(colHeader, collapse = ""))
    c(colHeader, fgs)
}



# /// Display functions (numeric GTs) ///////////////////////////////

## ----
# Format Reference Probability Values
#
# @description
# Formats reference probability values with color-coded backgrounds
#
# @param rpVals
# Numeric vector of reference probability values between 0 and 1
#
# @return
# Character vector of ANSI-formatted strings
formatRefProb <- function(rpVals) {
    if (!is.numeric(rpVals)) {
        rlang::abort("Input must be numeric.")
    }

    if (any(rpVals < 0 | rpVals > 1, na.rm = TRUE)) {
        rlang::abort("All values must be in the range [0–1]")
    }

    # ANSI background colors: white to turquoise (blue/green)
    bgCodes <- c(231, 195, 159, 87, 45)

    vapply(rpVals, FUN.VALUE = character(1), FUN = function(val) {
        if (is.na(val)) return(" NA ")

        idx <- ceiling(val / 0.2)
        idx <- max(1L, min(5L, idx)) # ensure 0 vals are 1 for indexing

        sprintf("\033[30;48;5;%dm %.3f \033[0m", bgCodes[idx], val)
    })
}


## ----
# Generate Numeric Genotype Display Strings
#
# @description
# Creates formatted display strings for numeric genotype data
#
# @param gt
# Genotype object
# @param nTaxa
# Maximum number of taxa to display
# @param nSites
# Maximum number of sites to display
#
# @return
# A list of formatted strings for display
genNumGtDispStrings <- function(gt, nTaxa = 5, nSites = 10) {
    totalSites <- gt$numberOfSites()
    totalTaxa  <- gt$numberOfTaxa()

    buildRowString <- function(taxonIndex, siteLimit, totalSites, gt) {
        if (totalSites > siteLimit + 2) {
            row <- vapply(seq_len(siteLimit), function(site) {
                formatRefProb(
                    gt$referenceProbability(
                        as.integer(taxonIndex - 1),
                        as.integer(site - 1)
                    )
                )
            }, FUN.VALUE = character(1))
            row <- c(
                row,
                subEllipsis(1),
                formatRefProb(
                    gt$referenceProbability(
                        as.integer(taxonIndex - 1),
                        as.integer(totalSites - 1)
                    )
                )
            )
        } else {
            row <- vapply(seq_len(totalSites), function(site) {
                formatRefProb(
                    gt$referenceProbability(
                        as.integer(taxonIndex - 1),
                        as.integer(site - 1)
                    )
                )
            }, FUN.VALUE = character(1))
        }
        return(row)
    }

    if (totalTaxa > nTaxa + 2) {
        returnStrings <- vector("list", nTaxa + 2)

        for (i in seq_len(nTaxa)) {
            returnStrings[[i]] <- buildRowString(i, nSites, totalSites, gt)
        }

        if (totalSites > nSites + 2) {
            ellipsisRow <- c(
                rep(subEllipsis(3), nSites),
                subEllipsis(1),
                subEllipsis(3)
            )
        } else {
            ellipsisRow <- rep(subEllipsis(3), totalSites)
        }

        returnStrings[[nTaxa + 1]] <- ellipsisRow
        returnStrings[[nTaxa + 2]] <- buildRowString(totalTaxa, nSites, totalSites, gt)

    } else {
        returnStrings <- vector("list", totalTaxa)
        for (i in seq_len(totalTaxa)) {
            returnStrings[[i]] <- buildRowString(i, nSites, totalSites, gt)
        }
    }

    return(returnStrings)
}


## ----
# Format Numeric Genotype Strings
#
# @description
# Formats numeric genotype data for display, including taxa names and
# IDs
#
# @param gt
# Genotype object
# @param nTaxa
# Maximum number of taxa to display
# @param nSites
# Maximum number of sites to display
#
# @return
# A list of formatted strings ready for display
formatNumGtStrings <- function(gt, nTaxa = 5, nSites = 10) {
    fGtStrings  <- genNumGtDispStrings(gt, nTaxa, nSites)

    fGtIds  <- addGtRowColIds(gt, fGtStrings, nSites, numeric = TRUE)
    fGtIds  <- lapply(fGtIds, paste, collapse = "")
    fGtTaxa <- formatGtNames(gt, nTaxa, 8)

    res <- Map(function(v, l) paste0(v, l), fGtTaxa, fGtIds)

    return(res)
}


## ----
# Print Numeric Genotype Display
#
# @description
# Prints formatted numeric genotype data to the console
#
# @param fgs
# Formatted genotype strings
# @param nTaxa
# Number of taxa
# @param nSites
# Number of sites
# @param jMem
# Java memory address
#
# @return
# None (called for side effects)
printNumGtDisp <- function(fgs, nTaxa, nSites, jMem) {
    header <- pillar::style_subtle(
        sprintf(
            "# A %s object: %s taxa %s %s sites\n",
            cli::style_bold("TasselNumericGenotype"),
            nTaxa,
            cli::symbol$times,
            nSites
        )
    )

    footer <- pillar::style_subtle(
        sprintf(
            "# %s Java memory address: 0x%s",
            cli::symbol$info,
            cli::style_bold(jMem)
        )
    )

    cat(header, "\n")
    for (i in seq_len(length(fgs))) {
        cat(fgs[[i]])
        cat("\n")
    }
    cat("\n")
    cat(footer, "\n")
}



# /// Display functions (allele GTs) ////////////////////////////////

## ----
# Format ANSI Bold Style
#
# @description
# Applies ANSI bold formatting to text
#
# @param allele
# Character string to format
#
# @return
# A character string with ANSI bold formatting
boldStyle <- function(allele) {
    sprintf("\033[1m %s \033[22m", allele)
}

## ----
# Format Green Background with Bold Text
#
# @description
# Applies green background with white bold text ANSI formatting
#
# @param allele
# Character string to format
#
# @return A character string with ANSI formatting
bgGreenBold <- function(allele) {
    sprintf("\033[42m\033[37m\033[1m %s \033[22m\033[39m\033[49m", allele)
}


## ----
# Format Yellow Background with Bold Text
#
# @description
# Applies yellow background with black bold text ANSI formatting
#
# @param allele
# Character string to format
#
# @return A character string with ANSI formatting
bgYellowBold <- function(allele) {
    sprintf("\033[43m\033[30m\033[1m %s \033[22m\033[39m\033[49m", allele)
}


## ----
# Format Blue Background with White Bold Text
#
# @description
# Applies blue background with white bold text ANSI formatting
#
# @param allele
# Character string to format
#
# @return A character string with ANSI formatting
bgBlueWhiteBold <- function(allele) {
    sprintf("\033[44m\033[37m\033[1m %s \033[22m\033[39m\033[49m", allele)
}


## ----
# Cache for allele Formatting Styles
styleCache <- list(
    "N"    = boldStyle("N"),
    "R"    = bgGreenBold("R"),
    "Y"    = bgGreenBold("Y"),
    "S"    = bgGreenBold("S"),
    "W"    = bgGreenBold("W"),
    "K"    = bgGreenBold("K"),
    "M"    = bgGreenBold("M"),
    "AMaj" = bgYellowBold("A"),
    "CMaj" = bgYellowBold("C"),
    "GMaj" = bgYellowBold("G"),
    "TMaj" = bgYellowBold("T"),
    "AMin" = bgBlueWhiteBold("A"),
    "CMin" = bgBlueWhiteBold("C"),
    "GMin" = bgBlueWhiteBold("G"),
    "TMin" = bgBlueWhiteBold("T")
)


## ----
# Format Allele with Styling
#
# @description
# This function formats a given allele based on its type and whether
# it matches the minimum allele. It applies specific styles from the
# `styleCache` object for recognized alleles and provides a default
# style for unrecognized ones.
#
# @details
# The function uses a `switch` statement to determine the appropriate
# style for the given allele. Recognized alleles include "A", "C",
# "G", "T", "N", "R", "Y", "S", "W", "K", "M", and "…". If the
# allele is not recognized, it is styled as bold with a default
# format.
#
# @param currAllele
# A character string representing the current allele.
# @param minAllele
# A character string representing the minimum allele for comparison.
#
# @return
# A styled character string corresponding to the formatted allele.
formatAllele <- function(currAllele, minAllele) {
    switch(currAllele,
        "A" = if (currAllele == minAllele) styleCache[["AMin"]] else styleCache[["AMaj"]],
        "C" = if (currAllele == minAllele) styleCache[["CMin"]] else styleCache[["CMaj"]],
        "G" = if (currAllele == minAllele) styleCache[["GMin"]] else styleCache[["GMaj"]],
        "T" = if (currAllele == minAllele) styleCache[["TMin"]] else styleCache[["TMaj"]],
        "N" = styleCache[["N"]],
        "R" = styleCache[["R"]],
        "Y" = styleCache[["Y"]],
        "S" = styleCache[["S"]],
        "W" = styleCache[["W"]],
        "K" = styleCache[["K"]],
        "M" = styleCache[["M"]],
        "…" = subEllipsis(1),
        sprintf("\033[1m %s \033[22m", currAllele) # default if not in cache
    )
}


## ----
# Generate Genotype Display Strings
#
# @description
# Generates formatted strings for displaying genotype data, handling
# truncation for large datasets
#
# @param gt
# Genotype object
# @param nTaxa
# Maximum number of taxa to display before truncating
# @param nSites
# Maximum number of sites to display before truncating
#
# @return
# A list of formatted genotype strings
genGtDispStrings <- function(gt, nTaxa = 5, nSites = 10) {
    maxTaxa <- gt$numberOfTaxa()
    maxSite <- gt$numberOfSites()

    getGt <- function(i) {
        i0 <- as.integer(i - 1)
        if (maxSite > nSites + 2) {
            head <- gt$genotypeAsStringRange(i0, 0L, as.integer(nSites))
            tail <- gt$genotypeAsString(i0, as.integer(maxSite - 1))
            paste0(head, cli::symbol$ellipsis, tail)
        } else {
            gt$genotypeAsStringRange(i0, 0L, as.integer(maxSite))
        }
    }

    if (maxTaxa > nTaxa + 2) {
        seqData <- c(
            seqData <- c(
                lapply(seq_len(nTaxa), getGt),
                list(
                    if (maxSite > nSites + 2) {
                        strrep(cli::symbol$ellipsis, nSites + 2)
                    } else {
                        strrep(cli::symbol$ellipsis, maxSite)
                    }
                ),
                list(getGt(maxTaxa))
            )
        )
    } else {
        seqData <- lapply(seq_len(maxTaxa), getGt)
    }

    lapply(seqData, function(it) unlist(strsplit(it, split = "")))
}


## ----
# Generate Minor Alleles List
#
# @description
# Creates a list of minor alleles for each site
#
# @param gt
# Genotype object
# @param nSites
# Maximum number of sites to process
#
# @return Vector of minor allele characters
genMinorAlleles <- function(gt, nSites = 10) {
    maxSite <- gt$numberOfSites()

    if (maxSite > nSites + 2) {
        minAlleles <- vapply(seq_len(nSites), function(it){
            gt$minorAlleleAsString(as.integer(it - 1))
        }, FUN.VALUE = character(1))
        minAlleles <- c(
            minAlleles,
            cli::symbol$ellipsis,
            gt$minorAlleleAsString(as.integer(gt$numberOfSites() - 1))
        )
    } else if (maxSite == nSites + 1) {
        minAlleles <- vapply(seq_len(nSites + 1), function(it){
            gt$minorAlleleAsString(as.integer(it - 1))
        }, FUN.VALUE = character(1))
    } else {
        minAlleles <- vapply(seq_len(maxSite), function(it){
            gt$minorAlleleAsString(as.integer(it - 1))
        }, FUN.VALUE = character(1))
    }

    return(minAlleles)
}


## ----
# Format Genotype Strings
#
# @description
# This function formats genotype strings for display by processing
# the input genotype data and applying various formatting operations.
# It generates display strings for genotypes, formats alleles based
# on minor alleles, and combines the results with row and column
# identifiers.
#
# @param gt
# A genotype object or data structure containing genotype
# information.
# @param nTaxa
# An integer specifying the number of taxa to include in the
# formatted output. Default is 5.
# @param nSites
# An integer specifying the number of sites to include in the
# formatted output. Default is 10.
#
# @details
# The function performs the following steps:
#   1. Generates display strings for genotypes using
#      `genGtDispStrings`.
#   2. Computes minor alleles for the given genotype data using
#      `genMinorAlleles`.
#   3. Formats each allele based on the minor allele information
#      using `formatAllele`.
#   4. Combines formatted genotype strings with row and column
#      identifiers using `addGtRowColIds`.
#   5. Formats taxa names for display using `formatGtNames`.
#   6. Combines formatted taxa names and genotype strings into the
#      final result.
#
# @return
# A list of formatted genotype strings, where each string represents
# a combination of formatted taxa names and genotype data.
formatGtStrings <- function(gt, nTaxa = 5, nSites = 10) {
    gtStrings  <- genGtDispStrings(gt, nTaxa, nSites)
    minAlleles <- genMinorAlleles(gt, nSites)

    fGtStrings <- lapply(gtStrings, function(gtString) {
        fgs <- Map(function(allele, ma) {
            formatAllele(allele, ma)
        }, gtString, minAlleles)
        paste0(fgs, collapse = "")
    })

    fGtIds  <- addGtRowColIds(gt, fGtStrings, nSites)
    fGtTaxa <- formatGtNames(gt, nTaxa, 8)

    res <- Map(function(v, l) paste0(v, l), fGtTaxa, fGtIds)

    return(res)
}


## ----
# Print Genotype Display
#
# @description
# This function prints a formatted display of genotype information,
# including the number of taxa, the number of sites, and a Java
# memory address. It also iterates through and prints the elements
# of the provided genotype data.
#
# @param fgs
# A list containing genotype data to be displayed.
# @param nTaxa
# An integer representing the number of taxa.
# @param nSites
# An integer representing the number of sites.
# @param jMem
# A string representing the Java memory address.
#
# @return
# This function does not return a value. It prints formatted genotype
# information to the console.
printGtDisp <- function(fgs, nTaxa, nSites, jMem) {
    header <- pillar::style_subtle(
        sprintf(
            "# A %s object: %s taxa %s %s sites\n",
            cli::style_bold("TasselGenotype"),
            nTaxa,
            cli::symbol$times,
            nSites
        )
    )

    footer <- pillar::style_subtle(
        sprintf(
            "# %s Java memory address: 0x%s",
            cli::symbol$info,
            cli::style_bold(jMem)
        )
    )

    cat(header, "\n")
    for (i in seq_len(length(fgs))) {
        cat(fgs[[i]])
        cat("\n")
    }
    cat("\n")
    cat(footer, "\n")
}



# /// Read functions ////////////////////////////////////////////////

## ----
#' Read Genotype Data from R Matrix
#'
#' @description
#' This function constructs a `TasselNumericGenotype` object from an
#' R matrix by interfacing with the TASSEL Java API. It creates taxa
#' lists, position lists, reference probabilities, and a genotype
#' table using the provided matrix data.
#'
#' @details
#' - Taxa list is created using the `TASSEL_JVM$TAXA_LIST_BUILDER`
#'   Java class.
#' - Position list is constructed using the
#'   `TASSEL_JVM$POSITION_LIST_BUILDER` and
#'   `TASSEL_JVM$GENERAL_POSITION_BUILDER` Java classes.
#' - Reference probabilities are built using the
#'   `TASSEL_JVM$REF_PROBABILITY_BUILDER` Java class.
#' - The genotype table is created using the
#'   `TASSEL_JVM$GENOTYPE_TABLE_BUILDER` Java class.
#'
#' @param m
#' A numeric matrix where rows represent taxa and columns represent
#' positions/sites. The matrix values are used to calculate reference
#' probabilities.
#' @param asTGP
#' Should the return object be a "classic" \code{TasselGenotypePhenotype}
#' object (\code{TRUE}) or should it return a \code{TasselNumericGenotype}
#' object (\code{FALSE})? Defaults to \code{TRUE}.
#'
#' @return
#' An object of class `TasselNumericGenotype` containing the genotype
#' data and metadata.
#'
#' @export
readNumericGenotypeFromRMatrix <- function(m, asTGP = TRUE) {
    if (!is(m, "matrix")) {
        rlang::abort("Provided object is not of type 'matrix'")
    }
    if (is.null(rownames(m))) {
        rlang::abort("No row names (taxa IDs) found")
    }
    if (is.null(colnames(m))) {
        rlang::abort("No column names (site IDs) found")
    }

    nRow <- nrow(m)
    nCol <- ncol(m)

    taxa <- rownames(m)
    mIds <- colnames(m)
    mPos <- seq_len(nCol)

    # Make taxa list
    tb <- rJava::.jnew(TASSEL_JVM$TAXA_LIST_BUILDER)
    taxaList <- tb$addAll(rJava::.jarray(taxa))$build()

    # Make position list
    pb <- rJava::.jnew(TASSEL_JVM$POSITION_LIST_BUILDER)
    chrUnknown <- rJava::.jnew(TASSEL_JVM$CHROMOSOME, "unknown")

    for (pos in mPos) {
        pb$add(
            rJava::.jnew(
                class = TASSEL_JVM$GENERAL_POSITION_BUILDER,
                chrUnknown,
                as.integer(pos - 1)
            )$snpName(
                mIds[pos]
            )$build()
        )
    }
    pbList <- pb$build()

    # Make ref. probabilities
    rpbClass <- rJava::J(TASSEL_JVM$REF_PROBABILITY_BUILDER)
    rpb <- rpbClass$getInstance(
        as.integer(nRow),
        as.integer(nCol),
        taxaList
    )
    for (i in seq_len(nRow)) {
        rpb$addTaxon(as.integer(i - 1), rJava::.jarray(rJava::.jfloat(m[i, ])))
    }
    rp <- rpb$build()

    # Make GenotypeTable
    gtbClass <- rJava::J(TASSEL_JVM$GENOTYPE_TABLE_BUILDER)
    javaGt <- gtbClass$getInstance(
        rJava::.jnull(), # genotype
        pbList,          # positions/sites
        taxaList,        # taxa
        rJava::.jnull(), # allele depth
        rJava::.jnull(), # allele probability
        rp,              # ref probability
        rJava::.jnull(), # dosage
        rJava::.jnull()  # annotations
    )

    jClass      <- rJava::.jclass(javaGt)
    jMemAddress <- gsub(".*@", "", rJava::.jstrVal(javaGt))

    if (asTGP) {
        gt <- .tasselObjectConstructor(javaGt)
    } else {
        gt <- methods::new(
            Class = "TasselNumericGenotype",
            dispData    = formatNumGtStrings(javaGt, nTaxa = 5, nSites = 5),
            jRefObj     = javaGt,
            jMemAddress = jMemAddress,
            jClass      = jClass
        )
    }

    return(gt)
}


## ----
# Format Genotype Strings
#
# @description
# This function formats genotype strings for display by processing
# the input genotype data.
#
# @details
# The function performs the following steps:
#   - Generates genotype display strings using `genGtDispStrings`.
#   - Computes minor alleles using `genMinorAlleles`.
#   - Formats alleles using `formatAllele` and combines them into a
#     single string.
#   - Adds row and column IDs to the formatted genotype strings using
#     `addGtRowColIds`.
#   - Formats the names of the taxa using `formatGtNames`.
#   - Combines the formatted taxa names and genotype strings into the
#     final result.
#
# @param gt
# A genotype object or data structure containing genotype
# information.
# @param nTaxa
# An integer specifying the number of taxa to include in the
# formatted output. Default is 5.
# @param nSites
# An integer specifying the number of sites to include in the
# formatted output. Default is 10.
#
# @return
# A list of formatted genotype strings, where each string represents
# a combination of taxa and site information.
readGenotypeFromPath <- function(x, sortPositions, keepDepth) {
    xNorm <- normalizePath(x, mustWork = FALSE)
    if (!file.exists(xNorm)) {
        rlang::abort("The input path is not a valid file")
    }

    rJc         <- rJava::.jnew(TASSEL_JVM$R_METHODS)
    javaGt      <- rJc$read(xNorm, keepDepth, sortPositions)
    jClass      <- rJava::.jclass(javaGt)
    jMemAddress <- gsub(".*@", "", rJava::.jstrVal(javaGt))

    if (javaGt$hasGenotype()) {
        methods::new(
            Class = "TasselGenotype",
            dispData    = formatGtStrings(javaGt),
            jRefObj     = javaGt,
            jMemAddress = jMemAddress,
            jClass      = jClass
        )
    } else {
        methods::new(
            Class = "TasselNumericGenotype",
            dispData    = formatNumGtStrings(javaGt, nTaxa = 5, nSites = 5),
            jRefObj     = javaGt,
            jMemAddress = jMemAddress,
            jClass      = jClass
        )
    }
}


