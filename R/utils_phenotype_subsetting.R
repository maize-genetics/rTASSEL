# /// Internal Helpers (phenotype metadata) //////////////////////////

## ----
# Wrapper classes that hold phenotype data in R alongside Java
#
# @description
# Both keep the tables that 'attributeData()' and 'as.data.frame()'
# return, so the metadata a verb needs can be read off the object rather
# than pulled back out of TASSEL.
PHENOTYPE_WRAPPERS <- c("TasselPhenotype", "TasselGenomicDataset")


## ----
# Attribute metadata for whichever phenotype a verb was handed
#
# @description
# The wrapper classes already carry the tibble that 'attributeData()'
# returns, so it is read off the object when possible and rebuilt from
# the Java 'Phenotype' only for the deprecated 'TasselGenotypePhenotype'.
#
# @param tasIn
# The list returned by '.resolveTasselInput()'.
#
# @return
# A 'tibble' with one row per attribute: 'trait_id', 'trait_type',
# 'trait_attribute', 'r_type', and the 0-based 'attr_idx'.
phenotypeAttrData <- function(tasIn) {
    if (.isAnyClass(tasIn$original, PHENOTYPE_WRAPPERS)) {
        return(attributeData(tasIn$original))
    }

    makeAttributeData(tasIn$jPh, tableReportToDF(tasIn$jPh))
}


## ----
# Observation data for whichever phenotype a verb was handed
#
# @param tasIn
# The list returned by '.resolveTasselInput()'.
#
# @return
# A 'tibble' of one row per observation, as 'as.data.frame()' returns.
phenotypeRowData <- function(tasIn) {
    if (methods::is(tasIn$original, "TasselPhenotype")) {
        return(tasIn$original@rData)
    }
    if (methods::is(tasIn$original, "TasselGenomicDataset")) {
        return(tasIn$original@phenotype@rData)
    }

    tableReportToDF(tasIn$jPh)
}


## ----
# Name of the column holding taxa IDs
#
# @description
# TASSEL names this attribute 'Taxa' when it reads a phenotype file, but
# a phenotype built from a 'data.frame' keeps the column name it was
# given, so it is always looked up by attribute type.
#
# @param attrData
# The tibble returned by 'phenotypeAttrData()'.
#
# @return
# A single 'character' value.
taxaColumnName <- function(attrData) {
    attrData$trait_id[attrData$trait_type == "taxa"][[1L]]
}


## ----
# Return all taxa names from a Java Phenotype
#
# @description
# The phenotype counterpart of 'batchTaxaNames()'. A phenotype may hold
# several observations of one taxon, so this is shorter than the taxa
# column of 'phenotypeRowData()' whenever taxa are replicated.
#
# @param jPh
# A Java 'Phenotype' object reference.
#
# @return
# A 'character' vector of taxa IDs.
phenotypeTaxaNames <- function(jPh) {
    rJava::J(TASSEL_JVM$R_METHODS)$
        genotypeTableToSampleNameArray(jPh$taxa())
}


## ----
# Attribute rows describing traits rather than taxa
#
# @description
# The taxa attribute is an axis, not a trait, so it is excluded from the
# trait verbs the way 'traitNames()' excludes it.
#
# @param attrData
# The tibble returned by 'phenotypeAttrData()'.
#
# @return
# A 'tibble', a subset of the rows of 'attrData'.
traitAttrRows <- function(attrData) {
    attrData[attrData$trait_type != "taxa", , drop = FALSE]
}


## ----
# Build the data mask a trait predicate is evaluated against
#
# @param traitRows
# The tibble returned by 'traitAttrRows()'.
# @param rData
# The tibble returned by 'phenotypeRowData()'.
#
# @return
# A 'tibble' with one row per trait, holding 'traitIndex' (1-based, over
# the traits only), 'traitId', 'traitType', and 'notMissing'.
buildTraitMetadata <- function(traitRows, rData) {
    tibble::tibble(
        traitIndex = seq_len(nrow(traitRows)),
        traitId    = traitRows$trait_id,
        traitType  = traitRows$trait_type,
        notMissing = vapply(
            traitRows$trait_id,
            function(id) mean(!is.na(rData[[id]])),
            FUN.VALUE  = numeric(1),
            USE.NAMES  = FALSE
        )
    )
}


## ----
# Columns the taxa data mask carries whatever class it was built from
TAXA_MASK_COLS <- c("taxaId", "notMissing", "het")


## ----
# Metadata columns that only a phenotype can supply
#
# @description
# 'filterTaxa()' on a 'TasselGenomicDataset' can be a genotype query or a
# phenotype query, and picks its path by looking for these names in the
# predicate. The shared mask columns are excluded because the genotype
# mask supplies them too, so a predicate naming only those keeps the
# genotype path and its plugin push-down.
#
# @param rData
# The tibble returned by 'phenotypeRowData()'.
#
# @return
# A 'character' vector of column names.
phenotypeOnlyVars <- function(rData) {
    setdiff(colnames(rData), TAXA_MASK_COLS)
}


## ----
# Build the data mask a phenotype taxa predicate is evaluated against
#
# @description
# The mask is the phenotype's own columns, one element per observation,
# plus a 'taxaId' alias for the taxa column so that a predicate written
# against a genotype table reads the same way here. On a
# 'TasselGenomicDataset' the genotype metrics are broadcast onto the
# observations of each taxon, so 'notMissing' and 'het' keep the meaning
# they have in 'taxaWhere()' and shadow any trait of the same name.
#
# @param tasIn
# The list returned by '.resolveTasselInput()'.
# @param rData
# The tibble returned by 'phenotypeRowData()'.
# @param vars
# Variable names the predicate refers to, from 'predicateVars()'. The
# genotype metrics cost one JNI round-trip per taxon, so they are only
# computed when asked for.
#
# @return
# A named 'list', suitable as the 'data' argument of 'evalPredicate()'.
phenotypeTaxaMask <- function(tasIn, rData, vars) {
    taxaCol <- taxaColumnName(phenotypeAttrData(tasIn))

    mask <- as.list(rData)
    mask$taxaId <- rData[[taxaCol]]

    needed <- intersect(vars, c("notMissing", "het"))
    if (length(needed) == 0L || rJava::is.jnull(tasIn$jGt)) {
        return(mask)
    }

    meta <- buildTaxaMetadata(tasIn$jGt, needed = needed)
    onto <- match(mask$taxaId, meta$taxaId)

    for (col in needed) {
        mask[[col]] <- meta[[col]][onto]
    }

    mask
}



# /// Internal Helpers (phenotype subsetting) ////////////////////////

## ----
# Keep a subset of a Java Phenotype's observations
#
# @description
# TASSEL has no row filter on 'PhenotypeBuilder', so each attribute is
# subset in turn and the results are reassembled. Passing the original
# type list along keeps every attribute's TASSEL type intact.
#
# @param jPh
# A Java 'Phenotype' object reference.
# @param positions
# A 1-based 'integer' vector of observations to keep.
#
# @return
# A Java 'Phenotype' object reference.
subsetPhenotypeObs <- function(jPh, positions) {
    jIdx   <- rJava::.jarray(as.integer(positions - 1L))
    jAttrs <- jPh$attributeListCopy()

    kept <- rJava::.jnew(TASSEL_JVM$ARRAY_LIST)
    for (i in seq_len(jAttrs$size()) - 1L) {
        jAttr <- rJava::.jcast(
            jAttrs$get(as.integer(i)),
            TASSEL_JVM$PHENO_ATTRIBUTE
        )
        kept$add(jAttr$subset(jIdx, jAttr$name()))
    }

    rJava::.jnew(TASSEL_JVM$PHENO_BUILDER)$
        fromAttributeList(kept, jPh$typeListCopy())$
        build()$
        get(0L)
}


## ----
# Keep a subset of a Java Phenotype's attributes
#
# @param jPh
# A Java 'Phenotype' object reference.
# @param attrIdx
# A 0-based 'integer' vector of attributes to keep, as held in the
# 'attr_idx' column of 'phenotypeAttrData()'. The taxa attribute is
# added if it is not already listed, since a phenotype without one
# cannot be joined to genotype data or read back as a data frame.
# @param taxaIdx
# The 0-based index of the taxa attribute.
#
# @return
# A Java 'Phenotype' object reference.
subsetPhenotypeTraits <- function(jPh, attrIdx, taxaIdx) {
    attrIdx <- sort(unique(c(as.integer(taxaIdx), as.integer(attrIdx))))

    rJava::.jnew(TASSEL_JVM$PHENO_BUILDER)$
        fromPhenotype(jPh)$
        keepAttributes(rJava::.jarray(attrIdx))$
        build()$
        get(0L)
}


## ----
# Does a taxa selector ask a question only a phenotype can answer?
#
# @description
# An object carrying both kinds of data has two row axes a taxa
# predicate could address, so the predicate decides: naming a phenotype
# column is a phenotype query, and anything else stays on the genotype
# table, where a simple threshold can be pushed down to a TASSEL filter
# plugin. Shared by 'filterTaxa()' and by '[' on a
# 'TasselGenomicDataset' so that the two grammars route alike.
#
# @param selector
# A 'TaxaSelector', or any other value a taxa index accepts.
# @param rData
# The tibble returned by 'phenotypeRowData()'.
#
# @return
# A single 'logical' value.
isPhenotypeTaxaPredicate <- function(selector, rData) {
    if (!methods::is(selector, "TaxaSelector")) return(FALSE)
    if (selector@type != "predicate") return(FALSE)

    any(predicateVars(selector@quo) %in% phenotypeOnlyVars(rData))
}


## ----
# Resolve a taxa selector to the taxa IDs it names
#
# @description
# The phenotype counterpart of the 'ids' branch of
# 'resolveTaxaIds()'. Predicates are handled by the caller, which has
# the data mask they are evaluated against.
#
# @param selector
# A 'character' vector or a 'TaxaSelector' of type '"ids"'.
#
# @return
# A 'character' vector of taxa IDs.
resolvePhenotypeTaxaIds <- function(selector) {
    if (is.character(selector)) return(selector)

    if (!methods::is(selector, "TaxaSelector")) {
        rlang::abort(
            "Taxa selector must be a character vector or TaxaSelector"
        )
    }

    if (selector@type == "ids") return(selector@ids)

    rlang::abort(paste0("Unknown TaxaSelector type: ", selector@type))
}


## ----
# Apply a taxa selector to a Java Phenotype
#
# @description
# The phenotype counterpart of 'applyTaxaSelector()'. Named taxa are
# kept whole, with every observation they have, while a predicate is
# evaluated one observation at a time, so a replicated taxon can pass
# on one row and fail on another.
#
# The positions are returned alongside the subset phenotype because a
# bracket call goes on to evaluate its trait index against the
# observations that survived this one.
#
# @param jPh
# A Java 'Phenotype' object reference.
# @param selector
# A 'character' vector or a 'TaxaSelector'.
# @param tasIn
# The list returned by '.resolveTasselInput()'.
# @param rData
# The tibble returned by 'phenotypeRowData()'.
#
# @return
# A 'list' with the subset Java 'Phenotype' 'jPh' and the 1-based
# 'integer' 'positions' of the observations kept.
applyPhenotypeTaxaSelector <- function(jPh, selector, tasIn, rData) {
    negate <- methods::is(selector, "TaxaSelector") && selector@negate

    isPredicate <- methods::is(selector, "TaxaSelector") &&
        selector@type == "predicate"

    positions <- if (isPredicate) {
        vars <- predicateVars(selector@quo)

        resolvePredicatePositions(
            selector@quo,
            phenotypeTaxaMask(tasIn, rData, vars),
            nrow(rData),
            "observations",
            negate = negate
        )
    } else {
        taxaCol <- taxaColumnName(phenotypeAttrData(tasIn))
        allTaxa <- unique(rData[[taxaCol]])
        ids     <- resolvePhenotypeTaxaIds(selector)

        if (negate) {
            ids <- setdiff(allTaxa, ids)
        } else {
            # An ID that names no taxon is skipped, as it is on a
            # genotype table, but naming nothing at all is an error
            if (length(setdiff(ids, allTaxa)) == length(ids)) {
                rlang::abort("No taxa match the selection criteria")
            }
            ids <- intersect(ids, allTaxa)
        }

        if (length(ids) == 0L) {
            rlang::abort("No taxa match the selection criteria")
        }

        which(rData[[taxaCol]] %in% ids)
    }

    list(
        jPh       = subsetPhenotypeObs(jPh, positions),
        positions = positions
    )
}


## ----
# Resolve a trait selector to 1-based positions among the traits
#
# @description
# Positions are counted as 'traitNames()' reports them, so the taxa
# column is not one of them.
#
# @param selector
# A 'numeric' vector, a 'character' vector, or a 'TraitSelector'.
# @param traitRows
# The tibble returned by 'traitAttrRows()'.
# @param rData
# The tibble returned by 'phenotypeRowData()'.
# @param negate
# Whether to keep the traits the selector did not name, as a selector
# inverted with '!' asks for.
#
# @return
# An 'integer' vector of 1-based positions in 'traitRows'.
resolveTraitPositions <- function(selector, traitRows, rData, negate = FALSE) {
    n <- nrow(traitRows)

    byName <- function(ids) {
        idx <- match(ids, traitRows$trait_id)
        as.integer(idx[!is.na(idx)])
    }

    byPosition <- function(idx) {
        idx <- as.integer(idx)

        if (any(idx < 1L) || any(idx > n)) {
            rlang::abort(c(
                "Trait positions must be 1-based and within the phenotype",
                "x" = sprintf("There %s %s trait%s to index",
                    if (n == 1L) "is" else "are", n, if (n == 1L) "" else "s"
                ),
                "i" = "`traitNames()` reports the traits in order"
            ))
        }

        idx
    }

    # A predicate takes its own complement, so that negating one that
    # matched nothing keeps every trait rather than raising an error
    complement <- function(idx) if (negate) setdiff(seq_len(n), idx) else idx

    if (is.numeric(selector)) return(complement(byPosition(selector)))
    if (is.character(selector)) return(complement(byName(selector)))

    if (!methods::is(selector, "TraitSelector")) {
        rlang::abort(
            "Trait selector must be numeric, character, or TraitSelector"
        )
    }

    switch(selector@type,
        "names" = complement(byName(selector@ids)),
        "predicate" = resolvePredicatePositions(
            selector@quo,
            buildTraitMetadata(traitRows, rData),
            n,
            "traits",
            negate = negate
        ),
        rlang::abort(paste0("Unknown TraitSelector type: ", selector@type))
    )
}


## ----
# Apply a trait selector to a Java Phenotype
#
# @description
# The trait-axis counterpart of 'applySiteSelector()'. The taxa column
# is the axis the observations sit on rather than a trait, so
# 'subsetPhenotypeTraits()' keeps it whatever the selector named.
#
# @param jPh
# A Java 'Phenotype' object reference.
# @param selector
# A 'numeric' vector, a 'character' vector, or a 'TraitSelector'.
# @param attrData
# The tibble returned by 'phenotypeAttrData()'.
# @param rData
# The tibble returned by 'phenotypeRowData()', holding the
# observations a predicate's 'notMissing' is computed over.
#
# @return
# A Java 'Phenotype' object reference.
applyTraitSelector <- function(jPh, selector, attrData, rData) {
    traitRows <- traitAttrRows(attrData)

    positions <- resolveTraitPositions(
        selector,
        traitRows,
        rData,
        negate = methods::is(selector, "TraitSelector") && selector@negate
    )

    if (length(positions) == 0L) {
        rlang::abort("No traits match the selection criteria")
    }

    subsetPhenotypeTraits(
        jPh,
        traitRows$attr_idx[positions],
        attrData$attr_idx[attrData$trait_type == "taxa"][[1L]]
    )
}


## ----
# Rebuild a phenotype result as the class that was handed in
#
# @description
# The phenotype counterpart of '.wrapGenotypeResult()'. Inputs that also
# carried genotype data are re-joined against the subset phenotype, which
# drops the taxa that no longer have any observations, so both components
# of the returned object stay in step.
#
# @param jPh
# The Java 'Phenotype' produced by a subsetting helper.
# @param tasIn
# The list returned by '.resolveTasselInput()' for the original input.
#
# @return
# An object of the same class as 'tasIn$original'.
.wrapPhenotypeResult <- function(jPh, tasIn) {
    jRes <- if (rJava::is.jnull(tasIn$jGt)) {
        jPh
    } else {
        joinGenotypePhenotype(tasIn$jGt, jPh)
    }

    .wrapLikeInput(jRes, tasIn$original)
}
