# /// Table report functions /////////////////////////////////////////

## ----
#' @title Table reports to tibble objects ----
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
    return(tibble::as_tibble(tabRepCols))
}


## ----
#' @title Get table reports based on HashMap ----
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


## ----
# Convert list to AssociationResults object ----
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



# /// LD plot functions //////////////////////////////////////////////

## ----
# Rotate vector coordinates by a given angle
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


## ----
# Rotated polygon coordinate, group, value "class"
ldCellRotater <- function(ldDF, angle) {
    # Reconstruct site rank order: coord2 first captures rank-1 site,
    # then coord1 adds the highest-ranked site not in coord2
    siteOrder <- unique(c(ldDF$coord2, ldDF$coord1))
    siteRank  <- stats::setNames(seq_along(siteOrder), siteOrder)

    n <- length(siteOrder) - 1L

    pairs <- which(lower.tri(matrix(0L, n, n), diag = TRUE), arr.ind = TRUE)
    pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]

    iIdx   <- pairs[, 1]
    jIdx   <- pairs[, 2]
    nCells <- nrow(pairs)
    nObs   <- nrow(ldDF)

    if (nObs == nCells) {
        fullVals <- ldDF[[3]]
    } else {
        # rank k (>= 2) → grid row k-1; rank j (>= 1) → grid col j
        obsRow <- siteRank[ldDF$coord1] - 1L
        obsCol <- siteRank[ldDF$coord2]
        obsPos <- obsRow * (obsRow - 1L) / 2L + obsCol

        fullVals <- rep(NA, nCells)
        fullVals[obsPos] <- ldDF[[3]]
    }

    # 4 vertices per cell (bl, tl, tr, br) in a single vectorized pass
    x     <- c(jIdx, jIdx, jIdx + 1, jIdx + 1)
    y     <- c(iIdx, iIdx + 1, iIdx + 1, iIdx)
    group <- rep.int(seq_len(nCells), 4L)
    val   <- rep.int(fullVals, 4L)

    rot <- rotate(x, y, angle)

    tibble::tibble(x = rot$x, y = rot$y, val = val, group = group)
}



# /// TasselDistanceMatrix functions /////////////////////////////////

## ----
# Truncate taxa IDs if too long
truncate <- function(t, max = 10, etc = "...") {
    if (nchar(t) > max) {
        return(paste0(c(unlist(strsplit(t, ""))[1:(max - 3)], etc), collapse = ""))
    } else {
        return(t)
    }
}


## ----
# Clean up taxa IDs and format spacing
cleanUpTaxa <- function(v, width = 10, regex = "^\"|\"$") {
    t <- gsub(regex, "", v)
    if (!is.null(width)) {
        t <- vapply(t, truncate, max = width, FUN.VALUE = "character")
        t <- format(t, width = width, justify = "right")
    }
    return(t)
}


## ----
# Clean up summary matrix with formatting
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


## ----
# "Pretty" print distance matrices
summaryDistance <- function(
    kinJ,
    width = 10,
    etc = "...",
    size = 5,
    regex = "^\"|\"$"
) {

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



# /// Data frame functions ///////////////////////////////////////////

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



# /// Class helper functions (TEMP) /////////////////////////////////

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



