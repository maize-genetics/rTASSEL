# /// S3 - helper functions /////////////////////////////////////////

## ----
# Create custom vector classes for phenotype data
covVctr  <- function(x, ...) vctrs::new_vctr(x, class = "cov")
dataVctr <- function(x, ...) vctrs::new_vctr(x, class = "data")
factVctr <- function(x, ...) vctrs::new_vctr(x, class = "fact")
taxaVctr <- function(x, ...) vctrs::new_vctr(x, class = "taxa")


## ----
# Create a Java-Compatible Phenotype Table
#
# @description
# This function creates a tibble from the provided data and assigns
# specific attributes to it, making it compatible with Java-based
# phenotype processing. The resulting object is assigned a custom
# class `"java_pheno_tbl"`.
#
# @param data
# A data frame or object that can be converted to a tibble.
# @param nTaxa
# An integer specifying the number of taxa.
# @param nTraits
# An integer specifying the number of traits.
# @param nCap
# An integer specifying the capacity (nCap) attribute.
# @param nDfRow
# An integer specifying the number of rows in the data frame.
# @param jMem
# A numeric value specifying the Java memory allocation (jMem).
#
# @return
# A tibble with additional attributes (`nCap`, `nTaxa`, `nTraits`,
# `nDfRow`, `jMem`) and a custom class `"java_pheno_tbl"`.
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
# Format Phenotype Display
#
# @description
# This function formats phenotype data for display by processing a
# data frame and its associated attribute metadata. It creates a
# table with formatted phenotype data and returns it as a
# Java-compatible object.
#
# @details
# The function processes the input data frame (`df`) and its
# associated metadata (`attrDf`) to create a formatted table. It
# uses the `trait_type` column in `attrDf` to determine how to
# process each trait in `df`. Supported trait types include:
#   - `"covariate"`...: Processed using `covVctr`.
#   - `"data"`........: Processed using `dataVctr`.
#   - `"factor"`......: Processed using `factVctr`.
#   - `"taxa"`........: Processed using `taxaVctr`.
#
# The resulting table is converted to a tibble and passed to
# `javaPhenoTbl` along with additional metadata such as the number
# of taxa, traits, and rows in the original data frame.
#
# @param df
# A data frame containing phenotype data.
# @param attrDf
# A data frame containing metadata about the traits in `df`. Each
# row should describe a trait with columns such as `trait_id` and
# `trait_type`.
# @param nCap
# An integer specifying the maximum number of rows to display.
# Defaults to 5.
# @param nTaxa
# An integer specifying the number of taxa in the dataset.
# @param jMem
# A Java memory object used for creating the Java-compatible table.
#
# @return
# A Java-compatible phenotype table object created using
# `javaPhenoTbl`.
formatPhenotypeDisplay <- function(df, attrDf, nCap = 5, nTaxa, jMem) {
    dfHead <- if (nrow(df) > nCap) head(df, nCap) else df

    tblData <- stats::setNames(
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
# Core trait selection method
#
# @description
# This function is the core component for downstream selection
# methods. It filters the attribute data based on the provided
# traits and the reference object. It ensures that the "Taxa" trait
# is always included in the selection. If any traits are missing,
# it will issue a warning. If no traits are found, it will abort
# the process with an error message.
#
# @param attrData
# A data frame or matrix containing attribute data.
# @param traits
# A character vector of trait names to be selected.
# @param jRefObj
# A reference object used to filter or compare traits.
#
# @return
# A subset of the attribute data containing only the selected common
# traits.
selectTraitsCommon <- function(attrData, traits, jRefObj) {
    # Ensure "Taxa" is in the trait list
    if (!"Taxa" %in% traits) {
        traits <- c("Taxa", traits)
    }

    # Filter attribute data for selected traits
    attrDataSub <- attrData[attrData$trait_id %in% traits, ]

    # Identify missing traits (excluding "Taxa" from the check)
    missingTraits <- setdiff(traits, attrData$trait_id)
    if (length(missingTraits) > 0) {
        rlang::warn(paste(
            "The following traits were not found in the phenotype data and will be ignored:",
            paste(missingTraits, collapse = ", ")
        ))
    }

    # Abort if no traits are found
    if (nrow(attrDataSub) == 0) {
        rlang::abort("No provided traits found in phenotype data")
    }

    # Build the phenotype using the Java builder
    phenoBuilder <- rJava::.jnew(TASSEL_JVM$PHENO_BUILDER)$
        fromPhenotype(jRefObj)$
        keepAttributes(rJava::.jarray(attrDataSub$attr_idx))$
        build()$
        get(0L)

    # Create and return the Tassel phenotype
    return(createTasselPhenotype(phenoBuilder))
}


## ----
# Select Traits from Formula
#
# @description
# This function selects specific traits from a phenotype object
# based on a given formula.
#
# @details
# The function first retrieves the attribute data from the phenotype
# object using `attributeData()`. It then parses the formula to
# identify the traits to keep using `parseFormula()`. Finally, it
# selects the specified traits from the phenotype object using
# `selectTraitsCommon()`.
#
# @param ph
# A phenotype object containing trait data.
# @param f
# A formula used to determine which traits to select.
#
# @return
# A subset of the phenotype object containing only the selected
# traits.
selectTraitsFromFormula <- function(ph, f) {
    attrDf <- attributeData(ph)

    traitsToKeep <- parseFormula(f, attrDf)

    return(
        selectTraitsCommon(attrDf, unlist(traitsToKeep), javaRefObj(ph))
    )
}



## ----
# Select Traits from Java Reference Object
#
# @description
# This function extracts and processes trait data from a Java
# reference object.
#
# @details
# The function first retrieves attribute data from the Java
# reference object by converting it into a data frame using
# `tableReportToDF` and `makeAttributeData`. It then uses a common
# helper function, `selectTraitsCommon`, to perform further
# processing and return the desired traits.
#
# @param jRefObj
# A Java reference object containing the data to be processed.
# @param traits
# A vector of trait names to be selected from the Java reference
# object.
#
# @return A processed data structure containing the selected traits.
selectTraitsFromJavaRef <- function(jRefObj, traits) {
    # Obtain necessary attribute data from the Java reference object
    attrData <- makeAttributeData(jRefObj, tableReportToDF(jRefObj))

    # Use the common helper for further processing
    return(selectTraitsCommon(attrData, traits, jRefObj))
}


## ----
# The TASSEL attribute types a phenotype column can be given
VALID_TASSEL_ATTRS <- c("taxa", "data", "covariate", "factor")


## ----
# Normalize an attribute type mapping
#
# @description
# Brings either accepted spelling of the `attrTypes` mapping to one
# canonical two-column tibble, which the rest of the data frame
# reader works from.
#
# @details
# Three forms are accepted:
#   - A named character vector, where the names are the columns of
#     the data frame and the values their TASSEL attribute types.
#     This is the hand-written form.
#   - A data frame with `col_id` and `tassel_attr` columns.
#   - A data frame with `trait_id` and `trait_type` columns, which is
#     what `attributeData()` reports, so the metadata read off one
#     phenotype can be used to build another.
#
# A character vector is only read by name, never by position, so a
# partially named or unnamed vector is rejected rather than silently
# matched up against the columns of the data frame.
#
# @param attrTypes
# A named character vector or a data frame in either spelling.
#
# @return
# A tibble with `col_id` and `tassel_attr` character columns.
normalizeAttrTypes <- function(attrTypes) {
    if (is.character(attrTypes)) {
        nms <- names(attrTypes)

        if (is.null(nms) || anyNA(nms) || !all(nzchar(nms))) {
            rlang::abort(c(
                "Every element of a character `attrTypes` needs a name",
                "i" = paste0(
                    "Name each element with the column it describes, e.g. ",
                    "c(Taxon = \"taxa\", EarHT = \"data\")."
                )
            ))
        }

        return(tibble::tibble(
            col_id      = nms,
            tassel_attr = unname(attrTypes)
        ))
    }

    if (is.data.frame(attrTypes)) {
        specCols <- c("col_id", "tassel_attr")
        attrCols <- c("trait_id", "trait_type")

        cols <- if (all(specCols %in% names(attrTypes))) {
            specCols
        } else if (all(attrCols %in% names(attrTypes))) {
            attrCols
        } else {
            rlang::abort(c(
                "`attrTypes` does not name its columns and attribute types",
                "i" = paste0(
                    "Use `col_id` and `tassel_attr`, or the `trait_id` and ",
                    "`trait_type` that `attributeData()` reports."
                )
            ))
        }

        return(tibble::tibble(
            col_id      = as.character(attrTypes[[cols[1]]]),
            tassel_attr = as.character(attrTypes[[cols[2]]])
        ))
    }

    rlang::abort(c(
        paste0(
            "`attrTypes` must be a named character vector or a `data.frame`, not <",
            class(attrTypes)[1], ">"
        ),
        "i" = paste0(
            "Either c(Taxon = \"taxa\", EarHT = \"data\") or a data frame with ",
            "`col_id` and `tassel_attr` columns."
        )
    ))
}


## ----
# Validate an attribute type mapping against its data frame
#
# @description
# Checks that a normalized mapping describes every column of the
# phenotype data frame exactly once, and that the TASSEL attribute
# types it names are legal.
#
# @details
# The mapping is held to the same contract `colData` is held to in a
# `SummarizedExperiment`: one entry per column, no entry without a
# column, and no column without an entry. Nothing is inferred from
# position or from the shape of the data.
#
# @param df
# The phenotype data frame.
# @param spec
# A mapping as returned by `normalizeAttrTypes()`.
#
# @return
# `spec`, invisibly. Called for the errors it raises.
validateAttrTypes <- function(df, spec) {
    dfCols <- names(df)

    if (nrow(df) == 0L) {
        rlang::abort("Phenotype `data.frame` has no observations")
    }

    if (is.null(dfCols) || anyNA(dfCols) || !all(nzchar(dfCols))) {
        rlang::abort("Every column of the phenotype `data.frame` needs a name")
    }

    dupCols <- unique(dfCols[duplicated(dfCols)])
    if (length(dupCols) > 0) {
        rlang::abort(paste0(
            "Phenotype `data.frame` has duplicate column names: ",
            paste(dupCols, collapse = ", ")
        ))
    }

    dupSpec <- unique(spec$col_id[duplicated(spec$col_id)])
    if (length(dupSpec) > 0) {
        rlang::abort(paste0(
            "`attrTypes` describes a column more than once: ",
            paste(dupSpec, collapse = ", ")
        ))
    }

    invalidAttrs <- setdiff(unique(spec$tassel_attr), VALID_TASSEL_ATTRS)
    if (length(invalidAttrs) > 0) {
        rlang::abort(c(
            paste0(
                "Illegal TASSEL attributes detected: ",
                paste(invalidAttrs, collapse = ", ")
            ),
            "i" = paste0(
                "Allowed attributes are: ",
                paste(VALID_TASSEL_ATTRS, collapse = ", ")
            )
        ))
    }

    nTaxaAttrs <- sum(spec$tassel_attr == "taxa")
    if (nTaxaAttrs != 1) {
        rlang::abort(c(
            paste0(
                "Exactly one 'taxa' attribute must be present in `attrTypes`, not ",
                nTaxaAttrs
            ),
            "i" = "The taxa column labels the observations, so there is always exactly one."
        ))
    }

    unmatched <- setdiff(spec$col_id, dfCols)
    undescribed <- setdiff(dfCols, spec$col_id)
    if (length(unmatched) > 0 || length(undescribed) > 0) {
        bullets <- character()
        if (length(unmatched) > 0) {
            bullets <- c(bullets, "x" = paste0(
                "Described in `attrTypes` but not a column of the data frame: ",
                paste(unmatched, collapse = ", ")
            ))
        }
        if (length(undescribed) > 0) {
            bullets <- c(bullets, "x" = paste0(
                "A column of the data frame but not described in `attrTypes`: ",
                paste(undescribed, collapse = ", ")
            ))
        }

        rlang::abort(c(
            "`attrTypes` must describe every column of the phenotype `data.frame`, and only those",
            bullets,
            "i" = "Add the column to `attrTypes`, or drop it from the data frame."
        ))
    }

    invisible(spec)
}


## ----
# Classify a phenotype column by its R storage
#
# @description
# Reports which of the R types a phenotype column can be built from
# a column is, or `NA` if it is none of them.
#
# @details
# Only bare atomic vectors and factors are accepted. Anything else -
# a `Date`, a `POSIXct`, a `difftime`, a list column, a matrix column
# - has no unambiguous TASSEL attribute, so it is rejected by name
# rather than being coerced on a guess.
#
# @param x
# A column of a phenotype data frame.
#
# @return
# One of "character", "double", "integer", "logical", or "factor",
# or `NA_character_`.
phenotypeColumnKind <- function(x) {
    if (is.factor(x)) return("factor")

    if (!is.null(oldClass(x)) || is.list(x) || !is.null(dim(x))) {
        return(NA_character_)
    }

    switch(typeof(x),
        "character" = "character",
        "double"    = "double",
        "integer"   = "integer",
        "logical"   = "logical",
        NA_character_
    )
}


## ----
# Coerce one phenotype column to what its TASSEL attribute needs
#
# @description
# Applies the one coercion rule that the column's TASSEL attribute
# type allows, and errors otherwise.
#
# @details
# TASSEL only builds a `NumericAttribute` from a `double[]` and a
# `CategoricalAttribute` from a `String[]`, so every column leaves
# here as one or the other:
#   - "taxa"......: character or factor, coerced with `as.character()`
#   - "data", "covariate": double or integer, coerced with `as.double()`
#   - "factor"....: anything but double, coerced with `as.character()`
#
# Missing values are legal in a numeric column, where TASSEL reads
# them as missing, but not in a taxa or factor column, where a taxon
# has to be named and `NA` would otherwise become a category of its
# own.
#
# @param x
# The column.
# @param colId
# The column's name, used in error messages.
# @param attrType
# The column's TASSEL attribute type.
#
# @return
# A `character` or `double` vector.
coercePhenotypeColumn <- function(x, colId, attrType) {
    kind <- phenotypeColumnKind(x)

    if (is.na(kind)) {
        rlang::abort(c(
            sprintf(
                "Column '%s' is of an unsupported type <%s>",
                colId, class(x)[1]
            ),
            "i" = paste0(
                "A phenotype column must be a character, double, integer, ",
                "logical, or factor vector."
            )
        ))
    }

    if (attrType == "taxa") {
        if (!kind %in% c("character", "factor")) {
            rlang::abort(c(
                sprintf(
                    "Column '%s' is marked as 'taxa' but is <%s>", colId, kind
                ),
                "i" = "Coerce it with `as.character()`."
            ))
        }

        out <- as.character(x)
        empty <- is.na(out) | !nzchar(trimws(out))
        if (any(empty)) {
            rlang::abort(c(
                sprintf(
                    "Column '%s' is marked as 'taxa' but has missing values", colId
                ),
                "i" = sprintf(
                    "First at row %s. Every observation needs a taxon ID.",
                    which(empty)[1]
                )
            ))
        }

        return(out)
    }

    if (attrType == "factor") {
        if (kind == "double") {
            rlang::abort(c(
                sprintf(
                    "Column '%s' is marked as 'factor' but is <double>", colId
                ),
                "i" = paste0(
                    "Coerce it with `as.character()`, or mark it as 'data' or ",
                    "'covariate'."
                )
            ))
        }

        out <- as.character(x)
        if (anyNA(out)) {
            rlang::abort(c(
                sprintf(
                    "Column '%s' is marked as 'factor' but has missing values", colId
                ),
                "i" = sprintf(
                    paste0(
                        "First at row %s. TASSEL reads every level of a factor as ",
                        "a category, so fill or drop these observations."
                    ),
                    which(is.na(out))[1]
                )
            ))
        }

        return(out)
    }

    if (!kind %in% c("double", "integer")) {
        rlang::abort(c(
            sprintf(
                "Column '%s' is marked as '%s' but is <%s>", colId, attrType, kind
            ),
            "i" = "Coerce it with `as.numeric()`, or mark it as 'factor'."
        ))
    }

    as.double(x)
}


## ----
# Coerce every phenotype column to what its TASSEL attribute needs
#
# @description
# Runs `coercePhenotypeColumn()` over the data frame, in the order
# the columns appear in it rather than the order `attrTypes` lists
# them.
#
# @param df
# The phenotype data frame.
# @param spec
# A mapping as returned by `normalizeAttrTypes()`, already checked
# against `df` by `validateAttrTypes()`.
#
# @return
# A named list of `character` and `double` vectors, in the column
# order of `df`.
coercePhenotypeColumns <- function(df, spec) {
    attrTypes <- stats::setNames(spec$tassel_attr, spec$col_id)

    stats::setNames(
        lapply(
            names(df),
            function(colId) {
                coercePhenotypeColumn(df[[colId]], colId, attrTypes[[colId]])
            }
        ),
        names(df)
    )
}


## ----
# Create Attribute Data Frame
#
# @description
# This function generates a data frame containing metadata about
# attributes (traits) from a Java object and corresponding R data.
# It merges information from both sources and returns a sorted
# tibble.
#
# @details
# The function performs the following steps:
#   - Extracts attribute metadata from the Java object using
#     \code{extractPhenotypeAttDf}.
#   - Appends R-side type information by applying the \code{class}
#     function to \code{rData}.
#   - Retrieves the attribute index from the Java object and creates
#     a cross-reference table.
#   - Merges the metadata and index data, sorts by attribute index,
#     and converts the result to a tibble.
#
# @param javaPh
# A Java object containing phenotype attribute information.
# It is expected to have methods for extracting attribute metadata
# and a list of attributes.
# @param rData
# An R object (e.g., a data frame or list) containing phenotype data.
# The function uses this to append R-side type information to the
# attribute metadata.
#
# @return
# A tibble containing the following columns:
#   - \code{trait_id}: The unique identifier for each trait.
#   - \code{trait_type}: The type of the trait (e.g., numeric, categorical).
#   - \code{trait_attribute}: Additional metadata about the trait.
#   - \code{r_type}: The R-side data type of the trait (e.g., "character", "numeric").
#   - \code{attr_idx}: The index of the trait as determined by the Java object.
# The tibble is sorted by the \code{attr_idx} column.
makeAttributeData <- function(javaPh, rData) {
    # Extract attribute metadata
    attrData <- extractPhenotypeAttDf(javaPh)

    # Append R-side type info
    attrData$r_type <- vapply(rData, class, "character")

    # Get trait index from Java side
    attrList <- rJava::.jevalArray(javaPh$attributeListCopy()$toArray())
    attrIdxXRef <- tibble::tibble(
        attr_idx = as.integer(seq_along(attrList) - 1),
        trait_id = .jStrings(attrList)
    )

    # Merge index data and return sorted df by attribute index
    attrData <- merge(attrData, attrIdxXRef, by = "trait_id")
    attrData <- attrData[order(attrData$attr_idx), ]

    tibble::as_tibble(attrData)
}


## ----
# Create a TasselPhenotype from a Java phenotype object
#
# @description
# This internal function builds a `TasselPhenotype` S4 object from a
# Java phenotype object.
#
# @param javaPh
# A Java object returned from the TASSEL phenotype builder.
#
# @return
# An S4 object of class `TasselPhenotype`.
createTasselPhenotype <- function(javaPh) {
    rData       <- tableReportToDF(javaPh)
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
# @description
# Internal helper function that reads and parses phenotype data from
# a file, constructs Java and R representations, and returns a
# `TasselPhenotype` object.
#
# @param
# path A character string representing the path to the phenotype file.
#
# @return
# An object of class `TasselPhenotype`.
readPhenotypeFromFile <- function(path) {
    xNorm <- normalizePath(path, mustWork = FALSE)
    if (!file.exists(xNorm)) {
        rlang::abort("The input path is not a valid file")
    }

    javaPh <- rJava::.jnew(TASSEL_JVM$PHENO_BUILDER)$fromFile(xNorm)$build()$get(0L)

    createTasselPhenotype(javaPh)
}


## ----
# Read Phenotype Data from Data Frame
#
# @description
# This function reads phenotype data from a data frame and a mapping
# of its columns to TASSEL attribute types, validates the input, and
# creates a TASSEL phenotype object.
#
# @details
# The function performs the following steps:
#   - Normalizes `attrTypes` to a canonical mapping.
#   - Checks that the mapping describes every column of `df` exactly
#     once, with a legal TASSEL attribute type.
#   - Coerces each column to the one R type its attribute allows, so
#     only `double[]` and `String[]` are handed to TASSEL.
#   - Constructs a TASSEL phenotype object from the results.
#
# @param df
# A data frame containing phenotype data. The data frame must include
# a column for taxa and other columns corresponding to phenotype
# attributes.
# @param attrTypes
# A named character vector, or a data frame in either spelling,
# mapping each column of `df` to its TASSEL attribute type. See
# `normalizeAttrTypes()`.
#
# @return
# A TASSEL phenotype object created from the input data frame.
readPhenotypeFromDf <- function(df, attrTypes) {
    spec <- normalizeAttrTypes(attrTypes)
    validateAttrTypes(df, spec)
    cols <- coercePhenotypeColumns(df, spec)

    specTypes <- stats::setNames(spec$tassel_attr, spec$col_id)
    taxaId    <- spec$col_id[spec$tassel_attr == "taxa"]
    traitIds  <- setdiff(names(cols), taxaId)

    # Make dataVector object of columns that are not "taxa" type
    dataVectors <- rJava::.jnew(TASSEL_JVM$ARRAY_LIST)
    for (traitId in traitIds) {
        dataVectors$add(rJava::.jarray(cols[[traitId]]))
    }

    # Make TASSEL phenotype object from valid dataframe components.
    # Each character vector is explicitly arrayed so that single-column
    # cases still resolve to the Java 'String[]' overload.
    rJc <- rJava::.jnew(TASSEL_JVM$R_METHODS)
    javaPh <- rJc$createPhenotypeFromRDataFrameElements(
        rJava::.jarray(cols[[taxaId]]),                      # taxa column
        rJava::.jarray(traitIds),                            # column IDs (not type taxa)
        rJava::.jarray(unname(specTypes[traitIds])),         # attribute types
        dataVectors                                          # non-taxa columns
    )

    createTasselPhenotype(javaPh)
}


