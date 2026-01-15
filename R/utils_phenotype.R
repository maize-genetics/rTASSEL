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
# `selectTraits()`.
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

    subPh <- selectTraits(ph, unlist(traitsToKeep))

    return(subPh)
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
# Select Traits from Phenotype Data
#
# @description
# This function retrieves attribute data from a phenotype object and
# processes it to select specific traits using a common helper
# function.
#
# @param ph
# A phenotype object containing the data to be processed.
# @param traits
# A vector of trait names to be selected from the phenotype data.
#
# @return
# A processed object containing the selected traits.
selectTraits <- function(ph, traits) {
    # Retrieve attribute data directly
    attrData <- attributeData(ph)

    # Use the common helper for further processing
    return(selectTraitsCommon(attrData, traits, javaRefObj(ph)))
}


## ----
# Validate Attribute Data Frame
#
# @description
# This function checks if the input `attrDf` is a valid data frame
# and contains the required columns: `col_id` and `tassel_attr`. If
# the input is not a data frame or if the required columns are
# missing, an error is raised.
#
# @param attrDf
# A data frame to be validated. It must contain the columns `col_id`
# and `tassel_attr`.
#
# @return
# This function does not return a value. It raises an error if the
# validation fails.
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
# Validate TASSEL Attributes
#
# @description
# This function validates the TASSEL attributes present in the
# provided data frames.
#
# @details
# The function performs the following validations:
#   - Ensures that all observed TASSEL attributes in `df` are valid.
#     The valid attributes are:
#       - "taxa"
#       - "covariate"
#       - "data"
#       - "factor"
#   - Ensures that the `attrDf` data frame contains exactly one
#     "taxa" attribute.
#
# If any invalid attributes are detected in `df`, or if `attrDf`
# does not contain exactly one "taxa" attribute, the function will
# throw an error using `rlang::abort`.
#
# @param df
# A data frame containing a column named `tassel_attr` with observed
# TASSEL attributes.
# @param attrDf
# A data frame containing a column named `tassel_attr` with
# attributes to validate.
#
# @return
# This function does not return a value. It is used for validation
# and will throw an error if the input data frames do not meet the
# specified criteria.
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
# Validate Columns in Data Frame
#
# @description
# This function checks whether the columns specified in the
# attribute data frame (`attrDf`) exist in the target data frame
# (`df`) and ensures that they are of the correct type.
#
# @details
# The function iterates through each row of `attrDf` and performs
# the following checks:
#   \itemize{
#     \item{Ensures that the column specified in `col_id` exists in `df`.}
#     \item{
#       Validates that columns marked as "data" or "covariate" in
#       `tassel_attr` are numeric.
#     }
#   }
# If any of these conditions are not met, the function throws an error
# using `rlang::abort`.
#
# @param df
# A data frame to be validated.
# @param attrDf
# A data frame containing column specifications. It must have
# two columns:
#   \describe{
#     \item{\code{col_id}}{The name of the column to validate in `df`.}
#     \item{\code{tassel_attr}}{The expected type of the column. It can be
#     either "data" or "covariate", which require the column to be numeric.}
#   }
#
# @return
# This function does not return a value. It is used for validation
# purposes and will throw an error if validation fails.
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

    buildResult <- rJava::.jnew(TASSEL_JVM$PHENO_BUILDER)$fromFile(xNorm)$build()
    javaPh <- safeGetFirst(buildResult)

    if (is.null(javaPh)) {
        abortPhenotypeLoadError(path)
    }

    createTasselPhenotype(javaPh)
}


## ----
# Read Phenotype Data from Data Frame
#
# @description
# This function reads phenotype data from a given data frame and an
# attribute data frame, validates the input, and creates a TASSEL
# phenotype object.
#
# @details
# The function performs the following steps:
#   - Validates the `attrDf` and cross-references it with `df` to
#     ensure consistency.
#   - Validates the TASSEL-specific attributes in `df` and `attrDf`.
#   - Creates a data vector object for columns that are not of type
#     "taxa".
#   - Constructs a TASSEL phenotype object using the validated data
#     frame components.
#
# @param df
# A data frame containing phenotype data. The data frame must include
# a column for taxa and other columns corresponding to phenotype
# attributes.
# @param attrDf
# A data frame describing the attributes of the columns in `df`. It
# must include the following columns:
#   - `col_id`: The column names in `df`.
#   - `tassel_attr`: The attribute type for each column (e.g., "taxa").
#
# @return
# A TASSEL phenotype object created from the input data frame.
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


