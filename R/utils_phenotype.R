# /// S3 - helper functions /////////////////////////////////////////

## ----
covVctr  <- function(x, ...) vctrs::new_vctr(x, class = "cov")
dataVctr <- function(x, ...) vctrs::new_vctr(x, class = "data")
factVctr <- function(x, ...) vctrs::new_vctr(x, class = "fact")
taxaVctr <- function(x, ...) vctrs::new_vctr(x, class = "taxa")


## ----
javaPhenoTbl <- function(data, nTaxa, nTraits, nCap, nDfRow, jMem) {
    df <- tibble::as_tibble(data)
    attr(df, "nCap")    <- nCap
    attr(df, "nTaxa")   <- nTaxa
    attr(df, "nTraits") <- nTraits
    attr(df, "nDfRow")  <- nDfRow
    attr(df, "jMem")    <- jMem

    class(df) <- c("java_pheno_tbl", class(df))

    df
}


## ----
formatPhenotypeDisplay <- function(df, attrDf, nCap = 5, nTaxa, jMem) {
    dfHead <- if (nrow(df) > nCap) head(df, nCap) else df

    tblData <- setNames(
        lapply(seq_len(nrow(attrDf)), function(i) {
            row <- attrDf[i, ]
            val <- dfHead[[row$trait_id]]
            switch(row$trait_type,
                "covariate" = covVctr(val),
                "data"      = dataVctr(val),
                "factor"    = factVctr(val),
                "taxa"      = taxaVctr(val)
            )
        }),
        attrDf$trait_id
    )

    return(
        javaPhenoTbl(
            data    = tibble::as_tibble(tblData),
            nTaxa   = nTaxa,
            nTraits = ncol(df),
            nCap    = nCap,
            nDfRow  = nrow(df),
            jMem    = jMem
        )
    )
}



# /// General utilities /////////////////////////////////////////////

## ----
selectTraitsFromFormula <- function(ph, f) {
    attrDf <- attributeData(ph)

    traitsToKeep <- parseFormula(f, attrDf)

    subPh <- selectTraits(ph, unlist(traitsToKeep))

    return(subPh)
}

## ----
selectTraits <- function(ph, traits) {
    # Check for "Taxa" id - this is needed!
    if (!"Taxa" %in% traits) {
        traits <- c("Taxa", traits)
    }

    # Filter attribute data for selected traits
    attrData <- attributeData(ph)
    attrDataSub <- attrData[attrData$trait_id %in% traits, ]

    # Identify missing traits (excluding "Taxa" from the check)
    missingTraits <- setdiff(traits, attrData$trait_id)
    if (length(missingTraits) > 0) {
        rlang::warn(paste(
            "The following traits were not found in the phenotype data and will be ignored:",
            paste(missingTraits, collapse = ", ")
        ))
    }

    if (nrow(attrDataSub) == 0) {
        rlang::abort("No provided traits found in phenotype data")
    }

    phenoBuilder <- rJava::.jnew(rTASSEL:::TASSEL_JVM$PHENO_BUILDER)$
        fromPhenotype(javaRefObj(ph))$
        keepAttributes(rJava::.jarray(attrDataSub$attr_idx))$
        build()$
        get(0L)

    return(rTASSEL:::createTasselPhenotype(phenoBuilder))
}


## ----
validateAttrDf <- function(attrDf) {
    # Ensure attrDf is a data frame and has required columns
    if (!inherits(attrDf, "data.frame")) {
        rlang::abort("'attrDf' parameter needs to be of type 'data.frame'")
    }
    requiredCols <- c("col_id", "tassel_attr")
    missingCols <- setdiff(requiredCols, names(attrDf))
    if (length(missingCols) > 0) {
        rlang::abort("Incorrect column IDs used - must be of type 'col_id' and 'tassel_attr'")
    }
}


## ----
validateTasselAttributes <- function(df, attrDf) {
    validTasselAttrs <- c("taxa", "covariate", "data", "factor")

    # Validate that observed tassel attributes in df are valid.
    obsTasselAttrs <- df[["tassel_attr"]]
    invalidAttrs <- setdiff(unique(obsTasselAttrs), validTasselAttrs)
    if (length(invalidAttrs) > 0) {
        rlang::abort(
            paste0(
                "Illegal TASSEL attributes detected: ",
                paste(invalidAttrs, collapse = ", "),
                "\n  Allowed attributes are:\n",
                paste("    *", validTasselAttrs, collapse = "\n")
            )
        )
    }

    # Ensure attrDf contains exactly one 'taxa' attribute
    if (sum(attrDf$tassel_attr == "taxa") != 1) {
        rlang::abort("Exactly one 'taxa' attribute must be present in 'attrDf'")
    }
}


## ----
validateColumns <- function(df, attrDf) {
    # Validate that each column specified in attrDf exists in df and is of the correct type
    for (i in seq_len(nrow(attrDf))) {
        colId <- attrDf$col_id[i]
        attrType <- attrDf$tassel_attr[i]

        if (!colId %in% colnames(df)) {
            rlang::abort(paste0("Column '", colId, "' listed in 'attrDf' not found in 'df'"))
        }

        colVal <- df[[colId]]
        if (attrType %in% c("data", "covariate") && !is.numeric(colVal)) {
            rlang::abort(paste0("Column '", colId, "' is marked as '", attrType, "' but is not numeric"))
        }
    }
}


## ----
makeAttributeData <- function(javaPh, rData) {
    # Extract attribute metadata
    attrData <- extractPhenotypeAttDf(javaPh)
    colnames(attrData) <- c("trait_id", "trait_type", "trait_attribute")

    # Append R-side type info
    attrData$r_type <- vapply(rData, class, "character")

    # Get trait index from Java side
    attrList <- rJava::.jevalArray(javaPh$attributeListCopy()$toArray())
    attrIdxXRef <- data.frame(
        attr_idx = as.integer(seq_along(attrList) - 1),
        trait_id = vapply(attrList, function(it) it$toString(), character(1))
    )

    # Merge index data and return sorted df by attribute index
    attrData <- merge(attrData, attrIdxXRef, by = "trait_id")
    attrData <- attrData[order(attrData$attr_idx), ]

    tibble::as_tibble(attrData)
}


## ----
# Create a TasselPhenotype from a Java phenotype object
#
# This internal function builds a `TasselPhenotype` S4 object from a Java
# phenotype object.
#
# @param javaPh A Java object returned from the TASSEL phenotype builder.
#
# @return An S4 object of class `TasselPhenotype`.
createTasselPhenotype <- function(javaPh) {
    rData       <- tibble::as_tibble(tableReportToDF(javaPh))
    jClass      <- rJava::.jclass(javaPh)
    jMemAddress <- gsub(".*@", "", rJava::.jstrVal(javaPh))

    attrData    <- makeAttributeData(javaPh, rData)
    attrSummary <- as.list(table(attrData$trait_type))

    dispData <- formatPhenotypeDisplay(
        df     = rData,
        attrDf = attrData,
        nCap   = 10,
        nTaxa  = javaPh$taxa()$numberOfTaxa(),
        jMem   = jMemAddress
    )

    methods::new(
        Class       = "TasselPhenotype",
        attrData    = attrData,
        attrSummary = attrSummary,
        dispData    = dispData,
        rData       = rData,
        jRefObj     = javaPh,
        jMemAddress = jMemAddress,
        jClass      = jClass
    )
}


## ----
# Read Phenotype Data from File
#
# Internal helper function that reads and parses phenotype data from a file,
# constructs Java and R representations, and returns a `TasselPhenotype` object.
#
# @param path A character string representing the path to the phenotype file.
#
# @return An object of class `TasselPhenotype`.
readPhenotypeFromFile <- function(path) {
    xNorm <- normalizePath(path)
    if (!file.exists(xNorm)) {
        rlang::abort("The input path is not a valid file")
    }

    javaPh <- rJava::.jnew(TASSEL_JVM$PHENO_BUILDER)$fromFile(xNorm)$build()$get(0L)

    createTasselPhenotype(javaPh)
}


## ----
readPhenotypeFromDf <- function(df, attrDf) {
    # Validate data and attribute cross references
    validateAttrDf(attrDf)
    validateTasselAttributes(df, attrDf)
    validateColumns(df, attrDf)

    # Make dataVector object of columns that are not "taxa" type
    dataVectors <- rJava::.jnew(TASSEL_JVM$ARRAY_LIST)
    nonTaxaRows <- attrDf[attrDf$tassel_attr != "taxa", ]
    for (i in seq_len(nrow(nonTaxaRows))) {
        colData <- df[[nonTaxaRows$col_id[i]]]
        dataVectors$add(rJava::.jarray(colData))
    }

    # Make TASSEL phenotype object from valid dataframe components
    rJc <- rJava::.jnew(TASSEL_JVM$R_METHODS)
    javaPh <- rJc$createPhenotypeFromRDataFrameElements(
        df[[attrDf[attrDf$tassel_attr == "taxa", ]$col_id]], # taxa column
        nonTaxaRows$col_id,                                  # column IDs (not type taxa)
        nonTaxaRows$tassel_attr,                             # attribute types
        dataVectors                                          # non-taxa columns
    )

    createTasselPhenotype(javaPh)
}


