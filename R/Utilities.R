#--------------------------------------------------------------------
# Script Name:   Utilities.R
# Description:   Utility functions for rTASSEL
# Author:        Brandon Monier
# Created:       2020-06-24 at 09:29:59
# Last Modified: 2020-06-24 at 09:37:10
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions for
#    utility and house-keeping purposes of main exported functions.
#--------------------------------------------------------------------

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
    return(S4Vectors::DataFrame(tabRepCols))
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


## Pretty print distance matrices ----
summaryDistance <- function(dmJ, m = 6) {
    etc <- "..."
    width <- 9
    simpleMat <- matrix(NA, m + 1, m + 1)

    elements <- c(1:(m - 1), dmJ$numberOfTaxa())
    taxa <- sapply(elements, function(i) dmJ$getTaxon(as.integer(i - 1))$toString())

    elements <- c(1, elements)
    taxa <- c("", taxa)
    mTW <- max(nchar(taxa))
    cutoff <- 8

    if (mTW > cutoff) {
        taxa[nchar(taxa) > cutoff] <- paste0(
            strtrim(taxa[nchar(taxa) > cutoff], cutoff - 3),
            "..."
        )
    }

    taxaR <- paste0("[", format(taxa, width = cutoff, justify = "right"), "]")
    taxaC <- paste0("[", format(taxa, width = cutoff, justify = "left"), "]")

    for (i in 1:(m + 1)) {
        for (j in 1:(m + 1)) {
            simpleMat[i, j] <- format(
                round(
                    x = dmJ$getDistance(
                        as.integer(elements[j] - 1),
                        as.integer(elements[i] - 1)
                    ),
                    digits = 5
                ),
                nsmall = 5,
                width = mTW
            )
        }
    }

    simpleMat[1, ] <- taxaR
    simpleMat[, 1] <- taxaC
    simpleMat <- format(simpleMat, width = cutoff, justify = "right")
    simpleMat[m, ] <- format(etc, justify = "centre", width = cutoff + 2)
    simpleMat[, m] <- " ... "
    simpleMat[1, 1] <- format("", width = cutoff + 2)

    return(simpleMat)
    # for (i in 1:(m + 1)) {
    #     cat(" ", simpleMat[i, ])
    #     cat("\n")
    # }
}

