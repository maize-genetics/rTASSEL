# /// Shared helpers ////////////////////////////////////////////////

## ----
# Collect Java Phenotype objects from a list of rTASSEL objects
#
# @description
# The join functions accept any mix of objects that carry phenotype
# data, plus 'PCAResults', whose principal components TASSEL also
# models as a phenotype. Each element is validated and unwrapped, and
# the class the joined result should be returned as is reported back.
#
# @param x
# A list (or vector) of rTASSEL objects.
# @param fn
# Name of the calling function, used in errors and warnings.
#
# @return
# A 'list' with elements 'jPhenotypes' (a Java 'ArrayList' of Java
# 'Phenotype' objects) and 'legacy' (a 'logical' that is 'TRUE' when
# every data-bearing input was a 'TasselGenotypePhenotype').
.collectPhenotypes <- function(x, fn) {
    if (length(x) == 0) {
        rlang::abort(
            sprintf("`x` must hold at least one object for `%s()`", fn)
        )
    }

    jPhenotypes <- rJava::.jnew(TASSEL_JVM$ARRAY_LIST)

    for (obj in x) {
        if (methods::is(obj, "PCAResults")) {
            jPhenotypes$add(obj@jObj)
        } else {
            jPhenotypes$add(.resolveTasselInput(obj, "phenotype", fn)$jPh)
        }
    }

    # 'PCAResults' is an analysis result rather than a data object, so it
    # does not get a vote on which class the join returns
    isLegacy <- vapply(x, .isAnyClass, logical(1), TASSEL_INPUT$LEGACY)
    isModern <- vapply(x, .isAnyClass, logical(1), TASSEL_INPUT$MODERN)

    list(
        jPhenotypes = jPhenotypes,
        legacy      = any(isLegacy) && !any(isModern)
    )
}


## ----
# Join a collection of Java Phenotype objects
#
# @param coll
# The list returned by '.collectPhenotypes()'.
# @param how
# The 'PhenotypeBuilder' join method to call: '"intersectJoin"',
# '"unionJoin"', or '"concatenate"'.
#
# @return
# A 'TasselPhenotype', or a 'TasselGenotypePhenotype' if every
# data-bearing input was one.
.joinPhenotypes <- function(coll, how) {
    builder <- rJava::.jnew(TASSEL_JVM$PHENO_BUILDER)$
        fromPhenotypeList(coll$jPhenotypes)

    builder <- switch(
        how,
        "intersectJoin" = builder$intersectJoin(),
        "unionJoin"     = builder$unionJoin(),
        "concatenate"   = builder$concatenate()
    )

    jPh <- builder$build()$get(0L)

    if (coll$legacy) {
        return(.tasselObjectConstructor(jPh))
    }

    return(createTasselPhenotype(jPh))
}



# /// Join methods //////////////////////////////////////////////////

## ----
#' @title Intersect join phenotype tables
#'
#' @description Intersect join multiple phenotype objects based on \code{Taxa}
#'    column.
#'
#' @param x A list of rTASSEL objects containing a phenotype. Accepted
#'    classes are \code{\linkS4class{TasselPhenotype}},
#'    \code{\linkS4class{TasselGenomicDataset}},
#'    \code{\linkS4class{PCAResults}}, and the deprecated
#'    \code{TasselGenotypePhenotype}.
#'
#' @return A \code{\linkS4class{TasselPhenotype}} object, or a
#'    \code{TasselGenotypePhenotype} if every data-bearing input was one.
#'
#' @importFrom rJava .jnew
#'
#' @export
intersectJoin <- function(x) {
    .joinPhenotypes(.collectPhenotypes(x, "intersectJoin"), "intersectJoin")
}


## ----
#' @title Union join phenotype tables
#'
#' @description Union join multiple phenotype objects based on \code{Taxa}
#'    column.
#'
#' @param x A list of rTASSEL objects containing a phenotype. Accepted
#'    classes are \code{\linkS4class{TasselPhenotype}},
#'    \code{\linkS4class{TasselGenomicDataset}},
#'    \code{\linkS4class{PCAResults}}, and the deprecated
#'    \code{TasselGenotypePhenotype}.
#'
#' @return A \code{\linkS4class{TasselPhenotype}} object, or a
#'    \code{TasselGenotypePhenotype} if every data-bearing input was one.
#'
#' @importFrom rJava .jnew
#'
#' @export
unionJoin <- function(x) {
    .joinPhenotypes(.collectPhenotypes(x, "unionJoin"), "unionJoin")
}


## ----
#' @title Concatenate phenotype tables
#'
#' @description Concatenate (e.g. bind rows) multiple phenotype objects based
#'    on \code{Taxa} column.
#'
#' @param x A list of rTASSEL objects containing a phenotype. Accepted
#'    classes are \code{\linkS4class{TasselPhenotype}},
#'    \code{\linkS4class{TasselGenomicDataset}},
#'    \code{\linkS4class{PCAResults}}, and the deprecated
#'    \code{TasselGenotypePhenotype}.
#'
#' @return A \code{\linkS4class{TasselPhenotype}} object, or a
#'    \code{TasselGenotypePhenotype} if every data-bearing input was one.
#'
#' @importFrom rJava .jnew
#'
#' @export
concatenate <- function(x) {
    .joinPhenotypes(.collectPhenotypes(x, "concatenate"), "concatenate")
}


##----
#' @title Merge genotype tables
#'
#' @description
#' Merges multiple genotype tables together by site information
#'
#' @return
#' A \code{\linkS4class{TasselGenotype}} object, or a
#' \code{TasselGenotypePhenotype} if every input was one.
#'
#' @name mergeGenotypeTables
#' @rdname mergeGenotypeTables
#'
#' @param tasObjL A list of objects containing genotype data. Accepted
#'    classes are \code{\linkS4class{TasselGenotype}},
#'    \code{\linkS4class{TasselGenomicDataset}}, and the deprecated
#'    \code{TasselGenotypePhenotype}.
#'
#' @export
mergeGenotypeTables <- function(tasObjL) {
    mergeGtClassPath <- "net/maizegenetics/analysis/data/MergeGenotypeTablesPlugin"
    gtClassPath      <- "net/maizegenetics/dna/snp/GenotypeTable"
    frameClassPath   <- "java/awt/Frame"

    if (!is(tasObjL, "list")) {
        rlang::abort("`tasObjL` must be a list")
    }
    if (length(tasObjL) == 0) {
        rlang::abort("`tasObjL` must hold at least one object")
    }

    jGts <- lapply(tasObjL, function(obj) {
        .resolveTasselInput(obj, "genotype", "mergeGenotypeTables")$jGt
    })

    gtArray <- rJava::.jarray(
        x = jGts,
        contents.class = gtClassPath
    )

    mergeGtPlugin <- rJava::new(
        rJava::J(mergeGtClassPath),
        rJava::.jnull(frameClassPath),
        FALSE
    )

    mergedGt <- mergeGtPlugin$mergeGenotypeTables(gtArray)

    if (all(vapply(tasObjL, .isAnyClass, logical(1), TASSEL_INPUT$LEGACY))) {
        return(.tasselObjectConstructor(mergedGt))
    }

    return(createTasselGenotype(mergedGt))
}
