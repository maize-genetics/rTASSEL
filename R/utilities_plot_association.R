## ----
#' @title Check for valid traits in association results
#' @param trait A \code{character} vector of trait IDs
#' @param assocRes An object of type \code{AssociationResults}
traitValidityChecker <- function(trait, assocRes) {
    if (!is.null(trait) && !all(trait %in% traitNames(assocRes))) {
        if (any(trait %in% traitNames(assocRes))) {
            misfitTraits <- paste(
                trait[!trait %in% traitNames(assocRes)],
                collapse = ", "
            )
            warning(
                "Some traits not found in results and will not be plotted:\n",
                "    ", misfitTraits
            )
        } else {
            stop("No traits specified are found in results")
        }
    }
}


## ----
#' @title Sanity checker for association results
#' @param x A given column to check
#' @param assocStats A \code{data.frame} object containing association stats
#' @param neededCols A \code{character} vector containing desired columns
assocStatsColumnChecker <- function(x, assocStats, neededCols) {
    if (!x %in% colnames(assocStats)) {
        stop(
            "'", x, "' column not found in stats dataframe - need at least: ",
            paste(neededCols, collapse = ",")
        )
    }
}


