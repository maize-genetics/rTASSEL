## ----
#' @title Intersect join phenotype tables
#'
#' @description Intersect join multiple phenotype objects based on \code{Taxa}
#'    column.
#'
#' @param x A vector of phenotype objects.
#'
#' @importFrom rJava .jnew
#'
#' @export
intersectJoin <- function(x) {
    phenoBuilder <- rJava::.jnew("net.maizegenetics.phenotype.PhenotypeBuilder")
    phenotypes   <- rJava::.jnew("java.util.ArrayList")

    classes <- vapply(x, class, "character")

    if (all(classes == "TasselGenotypePhenotype")) {
        capture <- lapply(x, function(i) phenotypes$add(i@jPhenotypeTable))

        intersectPhenotype <- .tasselObjectConstructor(
            phenoBuilder$
                fromPhenotypeList(phenotypes)$
                intersectJoin()$
                build()$
                get(0L)
        )
    } else {
        gp <- x[classes == "TasselGenotypePhenotype"]
        capture <- lapply(gp, function(i) phenotypes$add(i@jPhenotypeTable))

        lpca <- x[classes == "jobjRef"]
        capture <- lapply(lpca, function(i) phenotypes$add(i))

        intersectPhenotype <- rTASSEL:::.tasselObjectConstructor(
            phenoBuilder$
                fromPhenotypeList(phenotypes)$
                intersectJoin()$
                build()$
                get(0L)
        )
    }


    return(intersectPhenotype)
}


## ----
#' @title Union join phenotype tables
#'
#' @description Union join multiple phenotype objects based on \code{Taxa}
#'    column.
#'
#' @param x A vector of phenotype objects.
#'
#' @importFrom rJava .jnew
#'
#' @export
unionJoin <- function(x) {
    phenoBuilder <- rJava::.jnew("net.maizegenetics.phenotype.PhenotypeBuilder")
    phenotypes   <- rJava::.jnew("java.util.ArrayList")

    classes <- vapply(x, class, "character")

    if (all(classes == "TasselGenotypePhenotype")) {
        capture <- lapply(x, function(i) phenotypes$add(i@jPhenotypeTable))

        unionPhenotype <- .tasselObjectConstructor(
            phenoBuilder$
                fromPhenotypeList(phenotypes)$
                unionJoin()$
                build()$
                get(0L)
        )
    } else {
        gp <- x[classes == "TasselGenotypePhenotype"]
        capture <- lapply(gp, function(i) phenotypes$add(i@jPhenotypeTable))

        lpca <- x[classes == "jobjRef"]
        capture <- lapply(lpca, function(i) phenotypes$add(i))

        unionPhenotype <- rTASSEL:::.tasselObjectConstructor(
            phenoBuilder$
                fromPhenotypeList(phenotypes)$
                unionJoin()$
                build()$
                get(0L)
        )
    }

    return(unionPhenotype)
}


## ----
#' @title Concatenate phenotype tables
#'
#' @description Concatenate (e.g. bind rows) multiple phenotype objects based
#'    on \code{Taxa} column.
#'
#' @param x A vector of phenotype objects.
#'
#' @importFrom rJava .jnew
#'
#' @export
concatenate <- function(x) {
    phenoBuilder <- rJava::.jnew("net.maizegenetics.phenotype.PhenotypeBuilder")
    phenotypes   <- rJava::.jnew("java.util.ArrayList")

    capture <- lapply(x, function(i) phenotypes$add(i@jPhenotypeTable))

    concatPhenotype <- .tasselObjectConstructor(
        phenoBuilder$
            fromPhenotypeList(phenotypes)$
            concatenate()$
            build()$
            get(0L)
    )

    return(concatPhenotype)
}


