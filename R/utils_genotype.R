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
# Character vector of cli-styled strings
formatRefProb <- function(rpVals) {
    if (!is.numeric(rpVals)) {
        rlang::abort("Input must be numeric.")
    }

    if (any(rpVals < 0 | rpVals > 1, na.rm = TRUE)) {
        rlang::abort("All values must be in the range of 0 to 1")
    }

    bgColors <- c("#FFFFFF", "#D7FFFF", "#AFFFFF", "#5FFFFF", "#00D7FF")
    bgStyles <- lapply(bgColors, function(hex) {
        cli::combine_ansi_styles(
            cli::make_ansi_style(hex, bg = TRUE),
            cli::col_black
        )
    })

    vapply(rpVals, FUN.VALUE = character(1), FUN = function(val) {
        if (is.na(val)) return(" NA ")

        idx <- ceiling(val / 0.2)
        idx <- max(1L, min(5L, idx))

        bgStyles[[idx]](sprintf(" %.3f ", val))
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





# /// Display functions (allele GTs) ////////////////////////////////

## ----
# Allele styling helpers using cli
#
# These are called at display time (not cached at load time) so that
# cli can correctly detect terminal color capabilities.
styleBold <- function(allele) {
    cli::style_bold(sprintf(" %s ", allele))
}

styleHetero <- function(allele) {
    style <- cli::combine_ansi_styles(cli::bg_green, cli::col_white, cli::style_bold)
    style(sprintf(" %s ", allele))
}

styleMajor <- function(allele) {
    style <- cli::combine_ansi_styles(cli::bg_yellow, cli::col_black, cli::style_bold)
    style(sprintf(" %s ", allele))
}

styleMinor <- function(allele) {
    style <- cli::combine_ansi_styles(cli::bg_blue, cli::col_white, cli::style_bold)
    style(sprintf(" %s ", allele))
}


## ----
# Format Allele with Styling
#
# @description
# Formats a given allele based on its type and whether it matches the
# minor allele. Styling is computed at call time so cli can detect
# terminal capabilities correctly.
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
        "A" = if (currAllele == minAllele) styleMinor("A") else styleMajor("A"),
        "C" = if (currAllele == minAllele) styleMinor("C") else styleMajor("C"),
        "G" = if (currAllele == minAllele) styleMinor("G") else styleMajor("G"),
        "T" = if (currAllele == minAllele) styleMinor("T") else styleMajor("T"),
        "N" = styleBold("N"),
        "R" = styleHetero("R"),
        "Y" = styleHetero("Y"),
        "S" = styleHetero("S"),
        "W" = styleHetero("W"),
        "K" = styleHetero("K"),
        "M" = styleHetero("M"),
        "\u2026" = subEllipsis(1),
        styleBold(currAllele)
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

    siteGt <- function(i0, siteIdx) {
        gt$genotypeAsString(i0, as.integer(siteIdx))
    }

    getGt <- function(i) {
        i0 <- as.integer(i - 1)
        if (maxSite > nSites + 2) {
            headChars <- vapply(seq_len(nSites) - 1L, function(s) {
                siteGt(i0, s)
            }, FUN.VALUE = character(1))
            tailChar <- siteGt(i0, maxSite - 1L)
            paste0(c(headChars, cli::symbol$ellipsis, tailChar), collapse = "")
        } else {
            chars <- vapply(seq_len(maxSite) - 1L, function(s) {
                siteGt(i0, s)
            }, FUN.VALUE = character(1))
            paste0(chars, collapse = "")
        }
    }

    if (maxTaxa > nTaxa + 2) {
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
# @param fgs
# A list containing formatted genotype data to be displayed.
# @param nTaxa
# An integer representing the number of taxa.
# @param nSites
# An integer representing the number of sites.
# @param jMem
# A string representing the Java memory address.
# @param className
# Class name to display in the header.
#
# @return
# This function does not return a value. It prints formatted genotype
# information to the console.
printGtDisp <- function(fgs, nTaxa, nSites, jMem, className = "TasselGenotype") {
    header <- pillar::style_subtle(paste0(
        "# A ", cli::style_bold(className),
        " object: ", nTaxa, " taxa ", cli::symbol$times, " ", nSites, " sites"
    ))
    footer <- pillar::style_subtle(paste0(
        "# ", cli::symbol$info, " Java memory address: 0x", cli::style_bold(jMem)
    ))

    cli::cat_line(header)
    cli::cat_line()
    for (row in fgs) {
        cli::cat_line(row)
    }
    cli::cat_line()
    cli::cat_line(footer)
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
#' The following components are constructed using respective Java builder classes:
#' \itemize{
#'   \item Taxa list is created using the \code{TASSEL_JVM$TAXA_LIST_BUILDER}
#'         Java class.
#'   \item Position list is constructed using the
#'         \code{TASSEL_JVM$POSITION_LIST_BUILDER} and
#'         \code{TASSEL_JVM$GENERAL_POSITION_BUILDER} Java classes.
#'   \item Reference probabilities are built using the
#'         \code{TASSEL_JVM$REF_PROBABILITY_BUILDER} Java class.
#'   \item The genotype table is created using the
#'         \code{TASSEL_JVM$GENOTYPE_TABLE_BUILDER} Java class.
#' }

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
            jRefObj     = javaGt,
            jMemAddress = jMemAddress,
            jClass      = jClass
        )
    }

    return(gt)
}


## ----
# Read Genotype Data from File Path
#
# @description
# Reads genotype data from a file path using the TASSEL Java API and
# returns either a TasselGenotype or TasselNumericGenotype object
# depending on whether the data contains discrete genotype calls.
#
# @param x
# A normalized file path to the genotype data file.
# @param sortPositions
# Logical indicating whether to sort positions.
# @param keepDepth
# Logical indicating whether to retain depth information.
#
# @return
# A TasselGenotype or TasselNumericGenotype object.
readGenotypeFromPath <- function(x, sortPositions, keepDepth) {
    rJc         <- rJava::.jnew(TASSEL_JVM$R_METHODS)
    javaGt      <- rJc$read(x, keepDepth, sortPositions)
    jClass      <- rJava::.jclass(javaGt)
    jMemAddress <- gsub(".*@", "", rJava::.jstrVal(javaGt))

    if (javaGt$hasGenotype()) {
        methods::new(
            Class = "TasselGenotype",
            jRefObj     = javaGt,
            jMemAddress = jMemAddress,
            jClass      = jClass
        )
    } else {
        methods::new(
            Class = "TasselNumericGenotype",
            jRefObj     = javaGt,
            jMemAddress = jMemAddress,
            jClass      = jClass
        )
    }
}


