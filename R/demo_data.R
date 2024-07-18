## ----
#' @title Load Demonstration Data
#'
#' @description
#' This function loads demonstration data for phenotype, genotype,
#' or combined data types.
#'
#' @param type
#' A character vector indicating the type of data to load.
#' Possible values are "phenotype", "genotype", and "combine".
#'
#' @return A data object containing the requested type of data.
#'
#' @examples
#' # Load phenotype data
#' phenotype_data <- loadDemoData("phenotype")
#'
#' # Load genotype data
#' genotype_data <- loadDemoData("genotype")
#'
#' # Load combined data
#' combined_data <- loadDemoData("combine")
#'
#' @export
loadDemoData <- function(type = c("phenotype", "genotype", "combine")) {
    type <- rlang::arg_match(type)

    tasData <- list(
        "pheno" = system.file(
            "extdata", "mdp_traits.txt",
            package = "rTASSEL"
        ),
        "geno" = system.file(
            "extdata",
            "mdp_genotype.hmp.txt",
            package = "rTASSEL"
        )
    )

    tasObj <- switch(
        EXPR        = type,
        "phenotype" = readPhenotypeFromPath(tasData$pheno),
        "genotype"  = readGenotypeTableFromPath(tasData$geno),
        "combine"   = readGenotypePhenotype(
            tasData$geno,
            tasData$pheno
        )
    )

    return(tasObj)
}


