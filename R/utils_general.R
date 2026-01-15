# === Java build result helper functions ============================

## Safely get first element from a Java list build result ----
## Returns NULL if the list is empty, otherwise returns the first element
safeGetFirst <- function(javaList) {
    if (is.null(javaList) || rJava::is.jnull(javaList)) {
        return(NULL)
    }
    if (javaList$size() == 0L) {
        return(NULL)
    }
    return(javaList$get(0L))
}


## Check if a Java phenotype object has no taxa ----
phenotypeHasNoTaxa <- function(javaPhenotype) {
    if (is.null(javaPhenotype) || rJava::is.jnull(javaPhenotype)) {
        return(TRUE)
    }
    taxaList <- javaPhenotype$taxa()
    if (rJava::is.jnull(taxaList)) {
        return(TRUE)
    }
    return(taxaList$numberOfTaxa() == 0L)
}


## Extract sample taxa names from a Java TaxaList for error messages ----
#' @importFrom rJava is.jnull
extractTaxaSample <- function(javaTaxaList, maxSamples = 5) {
    if (is.null(javaTaxaList) || rJava::is.jnull(javaTaxaList)) {
        return(character(0))
    }
    nTaxa <- javaTaxaList$numberOfTaxa()
    if (nTaxa == 0) {
        return(character(0))
    }
    nSamples <- min(maxSamples, nTaxa)
    taxa <- vapply(
        seq_len(nSamples) - 1L,
        function(i) javaTaxaList$taxaName(as.integer(i)),
        FUN.VALUE = character(1)
    )
    return(taxa)
}


## Format taxa list for display in error messages ----
formatTaxaForMessage <- function(taxa, label, maxShow = 5) {
    if (length(taxa) == 0) {
        return(paste0(label, ": (none found)"))
    }
    taxaStr <- paste0("'", taxa[seq_len(min(length(taxa), maxShow))], "'", collapse = ", ")
    if (length(taxa) > maxShow) {
        taxaStr <- paste0(taxaStr, ", ...")
    }
    return(paste0(label, ": ", taxaStr))
}


# === Error message helper functions ================================

## Abort with phenotype file load error ----
abortPhenotypeLoadError <- function(path) {
    rlang::abort(c(
        paste0("Failed to load phenotype file: '", path, "'"),
        "x" = "The phenotype builder returned an empty result.",
        "i" = "This may indicate:",
        " " = "- The file format is not recognized",
        " " = "- The file is empty or contains no valid phenotype data",
        "i" = "Please check that the file is in a valid TASSEL phenotype format."
    ))
}


## Abort with no common taxa error (for genotype-phenotype intersection) ----
abortNoCommonTaxaGenoPheno <- function(genotypeTable, phenotype) {
    genoTaxa <- extractTaxaSample(genotypeTable$taxa())
    phenoTaxa <- extractTaxaSample(phenotype$taxa())

    genoMsg <- formatTaxaForMessage(genoTaxa, "Genotype taxa")
    phenoMsg <- formatTaxaForMessage(phenoTaxa, "Phenotype taxa")

    rlang::abort(c(
        "No common taxa found between genotype and phenotype data.",
        "x" = "The intersection of taxa IDs produced an empty result.",
        "i" = "This usually means the taxa/sample names in your genotype file do not match those in your phenotype file.",
        " " = "",
        "i" = "Sample taxa names from each dataset:",
        " " = genoMsg,
        " " = phenoMsg,
        " " = "",
        "i" = "Please ensure taxa names match exactly between files (case-sensitive)."
    ))
}


## Abort with no common taxa error (for phenotype intersection join) ----
abortNoCommonTaxaPhenoJoin <- function(taxaSamplesList) {
    taxaMsgs <- vapply(seq_along(taxaSamplesList), function(i) {
        formatTaxaForMessage(taxaSamplesList[[i]], paste0("Input ", i, " taxa"))
    }, FUN.VALUE = character(1))

    errMsgs <- c(
        "Intersect join produced no results.",
        "x" = "No common taxa were found across all input phenotype tables.",
        "i" = "The intersection requires taxa names to match exactly (case-sensitive).",
        " " = "",
        "i" = "Sample taxa names from each input:"
    )
    for (msg in taxaMsgs) {
        errMsgs <- c(errMsgs, " " = msg)
    }
    errMsgs <- c(errMsgs, " " = "", "i" = "Please ensure taxa names are consistent across all inputs.")

    rlang::abort(errMsgs)
}


## Collect taxa samples from a list of rTASSEL objects ----
collectTaxaSamplesFromObjects <- function(objList) {
    lapply(objList, function(obj) {
        if (inherits(obj, "TasselGenotypePhenotype") && !rJava::is.jnull(obj@jPhenotypeTable)) {
            extractTaxaSample(obj@jPhenotypeTable$taxa())
        } else if (inherits(obj, "PCAResults") && !rJava::is.jnull(obj@jObj)) {
            extractTaxaSample(obj@jObj$taxa())
        } else {
            character(0)
        }
    })
}


# === Table report functions ========================================

## Table reports to S4Vectors::DataFrame objects ----
#' @importFrom rJava .jevalArray
#' @importFrom rJava J
#' @importFrom S4Vectors DataFrame
tableReportToDF <- function(x) {
    rJC <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    tabRep <- rJC$tableReportToVectors(x)

    tabRepCols <- lapply(tabRep$dataVector, rJava::.jevalArray)

    tabRepCols <- do.call("data.frame", c(tabRepCols, stringsAsFactors = FALSE))
    colnames(tabRepCols) <- tabRep$columnNames
    colnames(tabRepCols) <- gsub(" ", "_", colnames(tabRepCols))
    # return(S4Vectors::DataFrame(tabRepCols, check.names = FALSE))
    return(tabRepCols)
}


## Get table reports based on HashMap ----
#' @importFrom rJava .jrcall
#' @importFrom rJava .jstrVal
tableReportList <- function(x) {

    hashVectors <- rJava::.jrcall(x, "keySet")
    hashVectors <- lapply(hashVectors, rJava::.jstrVal)

    myList <- lapply(hashVectors, function(i) x$get(i))
    myList <- lapply(myList, tableReportToDF)

    names(myList) <- hashVectors
    return(myList)
}


## Convert list to AssociationResults object ----
# @param trl A tableReportList object
# @param aType Association type
tableReportListToAssociationResults <- function(trl, aType) {
    result <- switch (aType,
        "BLUE" = {
            methods::new(
                Class = "AssociationResultsBLUE",
                results = trl,
                traits = trl$BLUE_ANOVA$Trait,
                assocType = aType
            )
        },
        "GLM" = {
            methods::new(
                Class = "AssociationResultsGLM",
                results = trl,
                traits = unique(trl$GLM_Stats$Trait),
                assocType = aType
            )
        },
        "MLM" = {
            methods::new(
                Class = "AssociationResultsMLM",
                results = trl,
                traits = unique(trl$MLM_Stats$Trait),
                assocType = aType
            )
        },
        "FastAssoc" = {
            methods::new(
                Class = "AssociationResultsFast",
                results = trl,
                traits = unique(trl$FastAssociation$Trait),
                assocType = aType
            )
        },
        "Stepwise" = {
            methods::new(
                Class = "AssociationResultsStepwise",
                results = trl,
                traits = unique(trl$ANOVA_report$Trait),
                assocType = aType
            )
        },
        "default" = {
            NULL
        }
    )

    if (is.null(result)) {
        stop("Association Type ('aType') not defined")
    } else {
        return(result)
    }
}



# === LD plot functions =============================================

## Rotate vector coordinates by a given angle ----
rotate <- function(x, y, angle = 135) {
    rad <- angle * (pi / 180)
    new_x <- x * cos(rad) - y * sin(rad)
    new_y <- y * cos(rad) + x * sin(rad)

    return(
        list(
            "x" = new_x,
            "y" = new_y
        )
    )
}


## Polygon coordinate, group, and value "class" ----
cell <- function(i, j, group, val, w = 1) {
    data.frame(
        ##        bl, tl,    tr,    br
        x     = c(j , j    , j + w, j + w),
        y     = c(i , i + w, i + w, i    ),
        group = group,
        val   = val
    )
}


## Rotated polygon coordinate, group, value "class" ----
ldCellRotater <- function(ldDF, angle) {
    rows <- cols <- length(unique(ldDF$coord1))
    grid_data <- matrix(data = data.frame(), nrow = rows, ncol = cols)

    # Generate path aesthetics - only half of matrix
    it <- 1
    for (i in seq_len(rows)) {
        for (j in seq_len(cols)) {
            if (i >= j) {
                sub <- ldDF[it, ]
                grid_data[[i, j]] <- cell(i, j, paste0(i, ":", j), sub[, 3])
                it <- it + 1
            }
        }
    }

    # Convert to data frame
    grid_data <- sapply(
        grid_data,
        "[",
        simplify = FALSE
    )
    grid_data <- do.call("rbind", grid_data)

    # Rotate and update coordinates
    rot <- rotate(grid_data$x, grid_data$y, angle)
    grid_data$x <- rot$x
    grid_data$y <- rot$y

    return(grid_data)
}



# === TasselDistanceMatrix functions ================================

## Truncate taxa IDs if too long ----
truncate <- function(t, max = 10, etc = "...") {
    if (nchar(t) > max) {
        return(paste0(c(unlist(strsplit(t, ""))[1:(max - 3)], etc), collapse = ""))
    } else {
        return(t)
    }
}


## Clean up taxa IDs and format spacing ----
cleanUpTaxa <- function(v, width = 10, regex = "^\"|\"$") {
    t <- gsub(regex, "", v)
    if (!is.null(width)) {
        t <- vapply(t, truncate, max = width, FUN.VALUE = "character")
        t <- format(t, width = width, justify = "right")
    }
    return(t)
}


## Clean up summary matrix with formatting ----
cleanUpMatrix <- function(t, d, space = "...", size = 5, width = 10, nTaxa) {
    if (nTaxa <= size) {
        m2 <- rbind(t, d)
        m2 <- cbind(c(format(" ", width = width), t), m2)
    } else {
        vec <- c(1, length(t) + 1)

        space <- format(space, width = width, justify = "right")

        m1 <- matrix(space, nrow = nrow(d) + length(vec), ncol = ncol(d))
        m1[-vec, ] <- d

        m2 <- matrix(space, nrow = nrow(m1), ncol = ncol(m1) + length(vec))
        m2[, -vec] <- m1
        m2[1, ] <- m2[, 1] <- c(
            space, t[1:(length(t) - 1)],
            space, t[length(t)]
        )
        m2[, nrow(m2) - 1] <- "   ..."
        m2[1, 1] <- format(" ", width = width)
    }


    return(m2)
}


## "Pretty" print distance matrices ----
summaryDistance <- function(kinJ,
                            width = 10,
                            etc = "...",
                            size = 5,
                            regex = "^\"|\"$") {

    if (kinJ$numberOfTaxa() <= size) {
        taxaCleaned <- cleanUpTaxa(
            v = kinJ$getTableColumnNames()[-1],
            width = width,
            regex = regex
        )

        distMat <- apply(
            X = rJava::.jevalArray(kinJ$getDistances(), simplify = TRUE),
            MARGIN = 1,
            format, digits = 4, nsmall = 4, justify = "right", width = 10
        )
    } else {
        indexes <- c(1:4, kinJ$numberOfTaxa())
        taxaCleaned <- cleanUpTaxa(
            v = kinJ$getTableColumnNames()[indexes + 1],
            width = width,
            regex = regex
        )

        distMat <- matrix(nrow = length(indexes), ncol = length(indexes))

        for (i in seq_along(indexes)) {
            for (j in seq_along(indexes)) {
                distMat[i, j] <- format(
                    sprintf(
                        kinJ$getDistance(
                            as.integer(indexes[i] - 1),
                            as.integer(indexes[j] - 1)
                        ),
                        fmt = "%#.4f"
                    ),
                    justify = "right",
                    width = 10
                )
            }
        }
    }

    return(
        cleanUpMatrix(
            t     = taxaCleaned,
            d     = distMat,
            space = etc,
            width = width,
            size  = size, nTaxa = kinJ$numberOfTaxa()
        )
    )
}



# === Data frame functions ==========================================

## ----
# @title Check for valid columns in a data frame object
# @param x A given column to check
# @param assocStats A \code{data.frame} object containing association stats
# @param neededCols A \code{character} vector containing desired columns
checkForValidColumns <- function(assocStats, neededCols) {
    for (col in neededCols) {
        if (!col %in% colnames(assocStats)) {
            stop(
                "'", col, "' column not found in stats dataframe - need at least: ",
                paste(neededCols, collapse = ",")
            )
        }
    }
}



# === Class helper functions (TEMP) =================================

## ----
# @title Get report elements
returnReportElements <- function(
    assocRes,
    reportName,
    defaultCatchAll = "ALL",
    defaultReportElement
) {
        if (!is.character(reportName) && !is.null(reportName)) {
            stop("'reportName' must be of type 'character'")
        }

        if (is.null(reportName)) {
            return(assocRes@results[[defaultReportElement]])
        }

        if (toupper(reportName) == defaultCatchAll) {
            return(assocRes@results)
        }

        if(reportName %in% reportNames(assocRes)) {
            return(assocRes@results[[reportName]])
        } else {
            stop("Report ID not found in object")
        }
}


