# /// Display functions (core) //////////////////////////////////////

## ----
subEllipsis <- function(ind = NULL) {
    spacer <- ifelse(!is.null(ind), strrep(" ", ind), "")
    sprintf(paste0(spacer, "%s", spacer), pillar::style_subtle(cli::symbol$ellipsis))
}


## ----
truncateGtName <- function(id, maxLength = 5) {
    if (nchar(id) > maxLength) {
        paste0(substr(id, 1, maxLength - 1), cli::symbol$ellipsis)
    } else {
        id
    }
}


## ----
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
formatRefProb <- function(rpVals) {
    if (!is.numeric(rpVals)) {
        rlang::abort("Input must be numeric.")
    }

    if (any(rpVals < 0 | rpVals > 1, na.rm = TRUE)) {
        rlang::abort("All values must be in the range [0–1]")
    }

    # ANSI background colors: white to turquoise
    bgCodes <- c(231, 195, 159, 87, 45)

    vapply(rpVals, FUN.VALUE = character(1), FUN = function(val) {
        if (is.na(val)) return(" NA ")

        idx <- ceiling(val / 0.2)
        idx <- max(1L, min(5L, idx)) # ensure 0 vals are 1 for indexing

        sprintf("\033[30;48;5;%dm %.3f \033[0m", bgCodes[idx], val)
    })
}


## ----
genNumGtDispStrings <- function(gt, nTaxa = 5, nSites = 10) {
    totalSites <- gt$numberOfSites()
    totalTaxa  <- gt$numberOfTaxa()

    buildRowString <- function(taxonIndex, siteLimit, totalSites, gt) {
        if (totalSites > siteLimit + 2) {
            row <- vapply(seq_len(siteLimit), function(site) {
                formatRefProb(
                    gt$referenceProbability(as.integer(taxonIndex - 1), as.integer(site - 1))
                )
            }, FUN.VALUE = character(1))
            row <- c(
                row,
                subEllipsis(1),
                formatRefProb(gt$referenceProbability(as.integer(taxonIndex - 1), as.integer(totalSites - 1)))
            )
        } else {
            row <- vapply(seq_len(totalSites), function(site) {
                formatRefProb(
                    gt$referenceProbability(as.integer(taxonIndex - 1), as.integer(site - 1))
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
formatNumGtStrings <- function(gt, nTaxa = 5, nSites = 10) {
    fGtStrings  <- genNumGtDispStrings(gt, nTaxa, nSites)

    fGtIds  <- addGtRowColIds(gt, fGtStrings, nSites, numeric = TRUE)
    fGtIds  <- lapply(fGtIds, paste, collapse = "")
    fGtTaxa <- formatGtNames(gt, nTaxa, 8)

    res <- Map(function(v, l) paste0(v, l), fGtTaxa, fGtIds)

    return(res)
}


## ----
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
# Mini functions for CLI formatting
# NOTE: hard-coding ANSI escapes here since calling CLI commands in a list at
#       run-time does not evaluate the full ANSI string?
boldStyle <- function(allele) {
    sprintf("\033[1m %s \033[22m", allele)
}

bgGreenBold <- function(allele) {
    sprintf("\033[42m\033[37m\033[1m %s \033[22m\033[39m\033[49m", allele)
}

bgYellowBold <- function(allele) {
    sprintf("\033[43m\033[30m\033[1m %s \033[22m\033[39m\033[49m", allele)
}

bgBlueWhiteBold <- function(allele) {
    sprintf("\033[44m\033[37m\033[1m %s \033[22m\033[39m\033[49m", allele)
}


## ----
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
           cli::style_bold(paste0(" ", currAllele, " "))  # Default for any other allele
    )
}


## ----
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



# /// Core builders /////////////////////////////////////////////////

## ----
readGenotypeFromRMatrix <- function(m) {
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

    methods::new(
        Class = "TasselNumericGenotype",
        dispData    = formatNumGtStrings(javaGt, nTaxa = 5, nSites = 5),
        jRefObj     = javaGt,
        jMemAddress = jMemAddress,
        jClass      = jClass
    )
}


## ----
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









