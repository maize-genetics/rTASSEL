## ----
# @title Check for valid traits in association results
# @param trait A \code{character} vector of trait IDs
# @param assocRes An object of type \code{AssociationResults}
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


