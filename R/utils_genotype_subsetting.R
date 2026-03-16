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
#
# @param needed Character vector of column names the caller actually
#   requires.  When non-NULL only the listed (expensive) columns are
#   computed, avoiding unnecessary JNI round-trips.  siteIndex, chrom,
#   and pos are always included because they come from a single batch
#   Java call.
buildSiteMetadata <- function(jGt, needed = NULL) {
    nSites    <- jGt$numberOfSites()
    positions <- jGt$positions()
    idx       <- seq_len(nSites) - 1L

    jRc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    plArrays <- jRc$genotypeTableToPositionListOfArrays(positions)

    meta <- data.frame(
        siteIndex = idx,
        chrom     = plArrays$chromosomes,
        pos       = as.numeric(plArrays$startPos),
        stringsAsFactors = FALSE
    )

    want <- function(col) is.null(needed) || col %in% needed

    # Try single-call batch Java helper (available when TASSEL provides
    # GenerateRCode.siteMetadataToArrays).  Falls back to per-site
    # rJava calls when the method is absent.
    batchArrays <- tryCatch(
        jRc$siteMetadataToArrays(jGt),
        error = function(e) NULL
    )

    if (!is.null(batchArrays)) {
        if (want("siteId"))      meta$siteId      <- batchArrays$siteNames
        if (want("maf"))         meta$maf          <- batchArrays$mafs
        if (want("alleleCount")) meta$alleleCount  <- as.integer(batchArrays$alleleCounts)
        if (want("het"))         meta$het          <- batchArrays$hets
        if (want("isIndel"))     meta$isIndel      <- batchArrays$indels
        if (want("isBiallelic")) meta$isBiallelic  <- batchArrays$biallelic
        return(meta)
    }

    if (want("siteId")) {
        meta$siteId <- batchSiteNames(jGt)
    }

    if (want("maf")) {
        meta$maf <- vapply(idx, function(i) {
            jGt$minorAlleleFrequency(as.integer(i))
        }, FUN.VALUE = numeric(1))
    }

    if (want("alleleCount")) {
        meta$alleleCount <- vapply(idx, function(i) {
            jGt$totalGametesNonMissingForSite(as.integer(i))
        }, FUN.VALUE = integer(1))
    }

    if (want("het")) {
        nTaxa <- jGt$numberOfTaxa()
        meta$het <- vapply(idx, function(i) {
            jGt$heterozygousCount(as.integer(i)) / nTaxa
        }, FUN.VALUE = numeric(1))
    }

    if (want("isIndel")) {
        meta$isIndel <- vapply(idx, function(i) {
            jGt$isIndel(as.integer(i))
        }, FUN.VALUE = logical(1))
    }

    if (want("isBiallelic")) {
        meta$isBiallelic <- vapply(idx, function(i) {
            alleles <- jGt$alleles(as.integer(i))
            length(alleles) <= 2L
        }, FUN.VALUE = logical(1))
    }

    meta
}


## ----
# Return all taxa names from a GenotypeTable via a single batch JNI
# call, avoiding per-element round-trips.
batchTaxaNames <- function(jGt) {
    jRc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    jRc$genotypeTableToSampleNameArray(jGt$taxa())
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
        taxaIds <- batchTaxaNames(jGt)

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

    allTaxa <- batchTaxaNames(jGt)

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
# Return all site names from a GenotypeTable.  Currently iterates via
# rJava; swap for a batch GenerateRCode helper once available.
batchSiteNames <- function(jGt) {
    nSites <- jGt$numberOfSites()
    vapply(seq_len(nSites) - 1L, function(i) {
        jGt$siteName(as.integer(i))
    }, FUN.VALUE = character(1))
}


## ----
# Resolve a site selector to 0-based integer site indices.
resolveSiteIndices <- function(jGt, selector) {
    if (is.numeric(selector) || is.integer(selector)) {
        return(as.integer(selector))
    }

    if (is.character(selector)) {
        allNames <- batchSiteNames(jGt)
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
            allNames <- batchSiteNames(jGt)
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
            needed <- all.vars(rlang::quo_get_expr(selector@quo))
            meta <- buildSiteMetadata(jGt, needed = needed)
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
# Try to translate a predicate SiteSelector into FilterSiteBuilderPlugin
# parameters.  Returns the filtered Java GenotypeTable on success, or
# NULL when the expression cannot be mapped to plugin parameters.
tryPluginShortCircuit <- function(jGt, selector) {
    if (!methods::is(selector, "SiteSelector") || selector@type != "predicate") {
        return(NULL)
    }
    if (selector@negate) return(NULL)

    params <- list(
        siteMinAlleleFreq    = 0.0,
        siteMaxAlleleFreq    = 1.0,
        minHeterozygous      = 0.0,
        maxHeterozygous      = 1.0,
        siteMinCount         = 0L,
        removeSitesWithIndels = FALSE
    )

    # Walk the AST and fill params. Returns TRUE if the full expression
    # was consumed, FALSE if any sub-expression is unsupported.
    consume <- function(expr) {
        if (!is.call(expr)) return(FALSE)

        op <- as.character(expr[[1]])

        if (op == "&" && length(expr) == 3L) {
            return(consume(expr[[2]]) && consume(expr[[3]]))
        }

        if (op %in% c(">=", "<=", ">", "<") && length(expr) == 3L) {
            lhs <- expr[[2]]
            rhs <- expr[[3]]
            if (!is.name(lhs) || !is.numeric(rhs)) return(FALSE)
            var <- as.character(lhs)
            val <- as.numeric(rhs)

            if (var == "maf") {
                if (op %in% c(">=", ">")) { params$siteMinAlleleFreq <<- val; return(TRUE) }
                if (op %in% c("<=", "<")) { params$siteMaxAlleleFreq <<- val; return(TRUE) }
            }
            if (var == "het") {
                if (op %in% c(">=", ">")) { params$minHeterozygous <<- val; return(TRUE) }
                if (op %in% c("<=", "<")) { params$maxHeterozygous <<- val; return(TRUE) }
            }
            if (var == "alleleCount") {
                if (op %in% c(">=", ">")) { params$siteMinCount <<- as.integer(val); return(TRUE) }
            }
            return(FALSE)
        }

        if (op == "!" && length(expr) == 2L && is.name(expr[[2]])) {
            if (as.character(expr[[2]]) == "isIndel") {
                params$removeSitesWithIndels <<- TRUE
                return(TRUE)
            }
        }

        FALSE
    }

    expr <- rlang::quo_get_expr(selector@quo)
    if (!consume(expr)) return(NULL)

    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin"),
        rJava::.jnull(),
        FALSE
    )
    plugin$setParameter("siteMinCount",         toString(params$siteMinCount))
    plugin$setParameter("siteMinAlleleFreq",    toString(params$siteMinAlleleFreq))
    plugin$setParameter("siteMaxAlleleFreq",    toString(params$siteMaxAlleleFreq))
    plugin$setParameter("minHeterozygous",      toString(params$minHeterozygous))
    plugin$setParameter("maxHeterozygous",      toString(params$maxHeterozygous))
    plugin$setParameter("removeSitesWithIndels", toString(params$removeSitesWithIndels))

    result <- tryCatch(plugin$runPlugin(jGt), error = function(e) NULL)
    if (!inherits(result, "jobjRef")) return(NULL)
    result
}


## ----
# Apply a site selector to a Java GenotypeTable. Returns a filtered
# Java GenotypeTable.
applySiteSelector <- function(jGt, selector) {
    shortCircuit <- tryPluginShortCircuit(jGt, selector)
    if (!is.null(shortCircuit)) return(shortCircuit)

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
