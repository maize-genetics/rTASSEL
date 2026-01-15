## ----
#' @title Intersect join phenotype tables
#'
#' @description Intersect join multiple phenotype objects based on \code{Taxa}
#'    column.
#'
#' @param x A list of rTASSEL objects containing a phenotype.
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
    } else {
        gp <- x[classes == "TasselGenotypePhenotype"]
        capture <- lapply(gp, function(i) phenotypes$add(i@jPhenotypeTable))

        lpca <- x[classes == "PCAResults"]
        capture <- lapply(lpca, function(i) phenotypes$add(i@jObj))
    }

    buildResult <- phenoBuilder$
        fromPhenotypeList(phenotypes)$
        intersectJoin()$
        build()

    result <- safeGetFirst(buildResult)

    # Check if result is NULL/empty OR if the resulting phenotype has no taxa
    if (is.null(result) || phenotypeHasNoTaxa(result)) {
        taxaSamples <- collectTaxaSamplesFromObjects(x)
        abortNoCommonTaxaPhenoJoin(taxaSamples)
    }

    intersectPhenotype <- .tasselObjectConstructor(result)
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
    } else {
        gp <- x[classes == "TasselGenotypePhenotype"]
        capture <- lapply(gp, function(i) phenotypes$add(i@jPhenotypeTable))

        lpca <- x[classes == "PCAResults"]
        capture <- lapply(lpca, function(i) phenotypes$add(i@jObj))
    }

    buildResult <- phenoBuilder$
        fromPhenotypeList(phenotypes)$
        unionJoin()$
        build()

    result <- safeGetFirst(buildResult)

    if (is.null(result)) {
        rlang::abort(c(
            "Union join produced no results.",
            "x" = "The phenotype builder returned an empty result.",
            "i" = "Please check that the input phenotype tables contain valid data."
        ))
    }

    unionPhenotype <- .tasselObjectConstructor(result)
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

    buildResult <- phenoBuilder$
        fromPhenotypeList(phenotypes)$
        concatenate()$
        build()

    result <- safeGetFirst(buildResult)

    if (is.null(result)) {
        rlang::abort(c(
            "Concatenation produced no results.",
            "x" = "The phenotype builder returned an empty result.",
            "i" = "Please check that the input phenotype tables contain valid data and have compatible columns for concatenation."
        ))
    }

    concatPhenotype <- .tasselObjectConstructor(result)
    return(concatPhenotype)
}


##----
#' @title Merge genotype tables
#'
#' @description
#' Merges multiple genotype tables together by site information
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @name mergeGenotypeTables
#' @rdname mergeGenotypeTables
#'
#' @param tasObjL A list of objects of class \code{TasselGenotypePenotype}.
#'
#' @export
mergeGenotypeTables <- function(tasObjL) {
    mergeGtClassPath <- "net/maizegenetics/analysis/data/MergeGenotypeTablesPlugin"
    gtClassPath      <- "net/maizegenetics/dna/snp/GenotypeTable"
    frameClassPath   <- "java/awt/Frame"

    if (!is(tasObjL, "list")) {
        stop("`tasObjL` must be a list")
    }

    if (!all(unlist(lapply(tasObjL, is, "TasselGenotypePhenotype")))) {
        stop("All elements in `tasObjL` must be of type TasselGenotypePhenotype")
    }

    gtArray <- rJava::.jarray(
        x = lapply(tasObjL, getGenotypeTable),
        contents.class = gtClassPath
    )

    mergeGtPlugin <- rJava::new(
        rJava::J(mergeGtClassPath),
        rJava::.jnull(frameClassPath),
        FALSE
    )

    mergedGt <- .tasselObjectConstructor(
        mergeGtPlugin$mergeGenotypeTables(gtArray)
    )

    return(mergedGt)
}


