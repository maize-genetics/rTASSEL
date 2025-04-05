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
           "…" = sprintf(" %s ", cli::symbol$ellipsis),
           cli::style_bold(paste0(" ", currAllele, " "))  # Default for any other allele
    )
}


## ----
genGtDispStrings <- function(gt, nTaxa = 5, nSites = 10) {
    seqData <- vector("list", length = nTaxa + 2)
    for (i in seq_len(nTaxa)) {
        head <- gt$genotypeAsStringRange(as.integer(i - 1), as.integer(0), as.integer(nSites))
        tail <- gt$genotypeAsString(as.integer(i - 1), as.integer(gt$numberOfSites() - 1))
        full <- paste0(head, cli::symbol$ellipsis, tail)

        seqData[[i]] <- full
    }

    seqData[[nTaxa + 1]] <- strrep(cli::symbol$ellipsis, nSites + 2)

    seqData[[nTaxa + 2]] <- paste0(
        gt$genotypeAsStringRange(as.integer(gt$numberOfTaxa() - 1), as.integer(0), as.integer(nSites)),
        cli::symbol$ellipsis,
        gt$genotypeAsString(as.integer(gt$numberOfTaxa() - 1), as.integer(gt$numberOfSites() - 1))
    )

    return(lapply(seqData, function(it) unlist(strsplit(it, split = ""))))
}


## ----
genMinorAlleles <- function(gt, nSites = 10) {
    minAlleles <- vapply(seq_len(nSites), function(it){
        gt$minorAlleleAsString(as.integer(it - 1))
    }, FUN.VALUE = character(1))

    minAlleles <- c(
        minAlleles,
        cli::symbol$ellipsis,
        gt$minorAlleleAsString(as.integer(gt$numberOfSites() - 1))
    )

    return(minAlleles)
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
    taxaIds <- vapply(seq_len(nTaxa), function(it) {
        gt$taxaName(as.integer(it - 1))
    }, FUN.VALUE = character(1))

    taxaIds <- c(
        "",
        taxaIds,
        paste0(strrep(" ", maxLength - 1), cli::symbol$ellipsis),
        gt$taxaName(as.integer(gt$numberOfTaxa() - 1))
    )

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
addGtRowColIds <- function(gt, fgs, nTaxa = 5, nSites = 10) {
    colHeader <- sprintf("%2d ", seq_len(nSites) - 1)
    colHeader <- c(
        colHeader,
        sprintf(" %s ", cli::symbol$ellipsis),
        sprintf(" %d", gt$numberOfSites() - 1)
    )

    colHeader <- pillar::style_subtle(paste(colHeader, collapse = ""))

    fgs <- c(colHeader, fgs)

    return(fgs)
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

    fGtIds  <- addGtRowColIds(gt, fGtStrings, nTaxa, nSites)
    fGtTaxa <- formatGtNames(gt, nTaxa, nSites)

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
    cat("\n", footer, "\n")
}


