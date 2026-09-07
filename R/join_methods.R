# /// Shared helpers ////////////////////////////////////////////////

## ----
# Flatten plain lists out of a variadic join call
#
# @description
# The join functions used to take a single list of objects, so
# 'c(ph1, ph2)' and 'list(ph1, ph2)' are still common call styles.
# Splicing bare lists out of '...' keeps them working alongside the
# variadic form. Only objects whose class is exactly 'list' are
# spliced, so a 'data.frame' still reaches the input validation and is
# reported as unsupported rather than being taken apart column by
# column.
#
# @param dots
# The 'list(...)' captured by a join function.
#
# @return
# A flat, unnamed 'list' of objects.
.spliceJoinInputs <- function(dots) {
    out <- list()

    for (i in seq_along(dots)) {
        el <- dots[[i]]

        if (identical(class(el), "list")) {
            out <- c(out, .spliceJoinInputs(el))
        } else {
            out[length(out) + 1L] <- list(el)
        }
    }

    return(unname(out))
}


## ----
# Collect Java Phenotype objects from a list of rTASSEL objects
#
# @description
# The join functions accept any mix of objects that carry phenotype
# data, plus 'PCAResults', whose principal components TASSEL also
# models as a phenotype. Objects carrying only genotype data are set
# aside so that '.joinPhenotypes()' can attach the joined phenotype to
# them. Each element is validated and unwrapped, and the class the
# joined result should be returned as is reported back.
#
# @param x
# A list (or vector) of rTASSEL objects.
# @param fn
# Name of the calling function, used in errors and warnings.
# @param allowGenotype
# Whether genotype-only input is meaningful for the calling function.
#
# @return
# A 'list' with elements 'jPhenotypes' (a Java 'ArrayList' of Java
# 'Phenotype' objects), 'jGts' (a 'list' of Java 'GenotypeTable'
# objects), and 'legacy' (a 'logical' that is 'TRUE' when every
# data-bearing input was a 'TasselGenotypePhenotype').
.collectPhenotypes <- function(x, fn, allowGenotype = TRUE) {
    if (length(x) == 0) {
        rlang::abort(
            sprintf("`%s()` needs at least one object to join", fn)
        )
    }

    jPhenotypes <- rJava::.jnew(TASSEL_JVM$ARRAY_LIST)
    jGts        <- list()
    nPh         <- 0L

    for (obj in x) {
        if (methods::is(obj, "PCAResults")) {
            jPhenotypes$add(obj@jObj)
            nPh <- nPh + 1L
            next
        }

        tasIn <- .resolveTasselInput(obj, "any", fn)

        if (!rJava::is.jnull(tasIn$jPh)) {
            jPhenotypes$add(tasIn$jPh)
            nPh <- nPh + 1L
        } else {
            jGts[length(jGts) + 1L] <- list(tasIn$jGt)
        }
    }

    if (nPh == 0L) {
        rlang::abort(c(
            sprintf(
                "`%s()` needs at least one object with phenotype data",
                fn
            ),
            "x" = "Only genotype data was found in the input"
        ))
    }

    if (!allowGenotype && length(jGts) > 0L) {
        rlang::abort(c(
            sprintf("`%s()` does not accept genotype data", fn),
            "i" = "Attach genotype data with `readGenomicDataset()` instead."
        ))
    }

    if (length(jGts) > 1L) {
        rlang::abort(c(
            sprintf("`%s()` accepts at most one genotype-only object", fn),
            "i" = "Combine genotype tables with `mergeGenotypeTables()` first."
        ))
    }

    # 'PCAResults' is an analysis result rather than a data object, so it
    # does not get a vote on which class the join returns
    isLegacy <- vapply(x, .isAnyClass, logical(1), TASSEL_INPUT$LEGACY)
    isModern <- vapply(x, .isAnyClass, logical(1), TASSEL_INPUT$MODERN)

    list(
        jPhenotypes = jPhenotypes,
        jGts        = jGts,
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
# A 'TasselPhenotype', a 'TasselGenomicDataset' when genotype data was
# collected, or a 'TasselGenotypePhenotype' if every data-bearing input
# was one.
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

    if (length(coll$jGts) == 0L) {
        if (coll$legacy) {
            return(.tasselObjectConstructor(jPh))
        }

        return(createTasselPhenotype(jPh))
    }

    jGp <- joinGenotypePhenotype(
        coll$jGts[[1L]],
        jPh,
        join = if (how == "unionJoin") "union" else "intersect"
    )

    if (coll$legacy) {
        return(.tasselObjectConstructor(jGp))
    }

    return(createTasselGenomicDataset(jGp))
}



# /// Join methods //////////////////////////////////////////////////

## ----
#' @title Intersect join phenotype tables
#'
#' @description Intersect join any number of phenotype objects based on the
#'    \code{Taxa} column. If one of the objects carries only genotype data,
#'    the joined phenotype is attached to it and a
#'    \code{\linkS4class{TasselGenomicDataset}} is returned.
#'
#' @param ... Any number of rTASSEL objects containing a phenotype. Accepted
#'    classes are \code{\linkS4class{TasselPhenotype}},
#'    \code{\linkS4class{TasselGenomicDataset}},
#'    \code{\linkS4class{PCAResults}}, and the deprecated
#'    \code{TasselGenotypePhenotype}. At most one object carrying only
#'    genotype data (\code{\linkS4class{TasselGenotype}}) may also be given.
#'    Lists of objects are flattened, so earlier
#'    \code{intersectJoin(c(ph1, ph2))} style calls keep working.
#'
#' @return A \code{\linkS4class{TasselPhenotype}} object, or a
#'    \code{\linkS4class{TasselGenomicDataset}} if genotype data was given.
#'    Returns a \code{TasselGenotypePhenotype} if every data-bearing input
#'    was one.
#'
#' @examples
#' \dontrun{
#' # Merge several phenotype tables of covariates and traits
#' intersectJoin(ph1Cov, ph2Traits, ph3MoreTraits)
#'
#' # Attach a genotype table at the same time
#' intersectJoin(gt, ph1Cov, ph2Traits, ph3MoreTraits)
#' }
#'
#' @importFrom rJava .jnew
#'
#' @export
intersectJoin <- function(...) {
    .joinPhenotypes(
        .collectPhenotypes(.spliceJoinInputs(list(...)), "intersectJoin"),
        "intersectJoin"
    )
}


## ----
#' @title Union join phenotype tables
#'
#' @description Union join any number of phenotype objects based on the
#'    \code{Taxa} column. If one of the objects carries only genotype data,
#'    the joined phenotype is attached to it and a
#'    \code{\linkS4class{TasselGenomicDataset}} is returned.
#'
#' @param ... Any number of rTASSEL objects containing a phenotype. Accepted
#'    classes are \code{\linkS4class{TasselPhenotype}},
#'    \code{\linkS4class{TasselGenomicDataset}},
#'    \code{\linkS4class{PCAResults}}, and the deprecated
#'    \code{TasselGenotypePhenotype}. At most one object carrying only
#'    genotype data (\code{\linkS4class{TasselGenotype}}) may also be given.
#'    Lists of objects are flattened, so earlier
#'    \code{unionJoin(c(ph1, ph2))} style calls keep working.
#'
#' @return A \code{\linkS4class{TasselPhenotype}} object, or a
#'    \code{\linkS4class{TasselGenomicDataset}} if genotype data was given.
#'    Returns a \code{TasselGenotypePhenotype} if every data-bearing input
#'    was one.
#'
#' @examples
#' \dontrun{
#' # Merge several phenotype tables of covariates and traits
#' unionJoin(ph1Cov, ph2Traits, ph3MoreTraits)
#'
#' # Attach a genotype table at the same time
#' unionJoin(gt, ph1Cov, ph2Traits, ph3MoreTraits)
#' }
#'
#' @importFrom rJava .jnew
#'
#' @export
unionJoin <- function(...) {
    .joinPhenotypes(
        .collectPhenotypes(.spliceJoinInputs(list(...)), "unionJoin"),
        "unionJoin"
    )
}


## ----
#' @title Concatenate phenotype tables
#'
#' @description Concatenate (e.g. bind rows) any number of phenotype objects
#'    based on the \code{Taxa} column.
#'
#' @param ... Any number of rTASSEL objects containing a phenotype. Accepted
#'    classes are \code{\linkS4class{TasselPhenotype}},
#'    \code{\linkS4class{TasselGenomicDataset}},
#'    \code{\linkS4class{PCAResults}}, and the deprecated
#'    \code{TasselGenotypePhenotype}. Lists of objects are flattened, so
#'    earlier \code{concatenate(c(ph1, ph2))} style calls keep working.
#'    Unlike the joins, this function binds phenotype rows together and so
#'    does not accept genotype-only input.
#'
#' @return A \code{\linkS4class{TasselPhenotype}} object, or a
#'    \code{TasselGenotypePhenotype} if every data-bearing input was one.
#'
#' @examples
#' \dontrun{
#' # Stack phenotype tables that share the same traits
#' concatenate(ph2021, ph2022, ph2023)
#' }
#'
#' @importFrom rJava .jnew
#'
#' @export
concatenate <- function(...) {
    .joinPhenotypes(
        .collectPhenotypes(
            .spliceJoinInputs(list(...)),
            "concatenate",
            allowGenotype = FALSE
        ),
        "concatenate"
    )
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
