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
