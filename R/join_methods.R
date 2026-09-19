# /// Constants /////////////////////////////////////////////////////

## ----
# Analysis result classes that keep a TASSEL 'Phenotype' in '@jObj'
#
# @description
# These are not data objects, but the values they report - PCA and MDS
# axes, and BLUE trait estimates - are held as a TASSEL phenotype, so
# the join functions can treat them as another table of traits.
PHENOTYPE_BACKED_RESULTS <- c(
    "PCAResults",
    "MDSResults",
    "AssociationResultsBLUE"
)



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
# Collect the Java data objects out of a list of rTASSEL objects
#
# @description
# The join functions accept any mix of objects that carry phenotype
# data, plus the analysis results whose values TASSEL also models as a
# phenotype: the axes of 'PCAResults' and 'MDSResults' and the trait
# estimates of 'AssociationResultsBLUE'. Objects carrying only genotype
# data are set aside so that '.joinPhenotypes()' can combine them and
# attach any joined phenotype to the result. Each element is validated
# and unwrapped, and the class the joined result should be returned as
# is reported back.
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

    for (obj in x) {
        if (.isAnyClass(obj, PHENOTYPE_BACKED_RESULTS)) {
            if (rJava::is.jnull(obj@jObj)) {
                rlang::abort(c(
                    sprintf(
                        "`%s()` got a <%s> object with no phenotype data",
                        fn, class(obj)
                    ),
                    "i" = paste0(
                        "Only results returned by the analysis functions ",
                        "carry a TASSEL phenotype."
                    )
                ))
            }

            jPhenotypes$add(obj@jObj)
            next
        }

        if (methods::is(obj, "AssociationResults")) {
            rlang::abort(c(
                sprintf(
                    "`%s()` cannot join <%s> results", fn, class(obj)
                ),
                "x" = paste0(
                    "Only BLUE results hold values TASSEL models as a ",
                    "phenotype - the other models report marker statistics."
                )
            ))
        }

        tasIn <- .resolveTasselInput(obj, "any", fn)

        if (!rJava::is.jnull(tasIn$jPh)) {
            jPhenotypes$add(tasIn$jPh)
        } else {
            jGts[length(jGts) + 1L] <- list(tasIn$jGt)
        }
    }

    if (!allowGenotype && length(jGts) > 0L) {
        rlang::abort(c(
            sprintf("`%s()` does not accept genotype data", fn),
            "i" = "Attach genotype data with `readGenomicDataset()` instead."
        ))
    }

    # The classes in 'PHENOTYPE_BACKED_RESULTS' are analysis results rather
    # than data objects, so they do not get a vote on the returned class
    isLegacy <- vapply(x, .isAnyClass, logical(1), TASSEL_INPUT$LEGACY)
    isModern <- vapply(x, .isAnyClass, logical(1), TASSEL_INPUT$MODERN)

    list(
        jPhenotypes = jPhenotypes,
        jGts        = jGts,
        legacy      = any(isLegacy) && !any(isModern)
    )
}


## ----
# Combine a collection of Java GenotypeTable objects
#
# @description
# Sites are taken from every table and taxa are combined the way the
# calling join asks for: an intersect keeps the taxa every table holds,
# a union keeps the taxa any of them holds and fills the calls the
# others never made with missing.
#
# TASSEL's 'CombineGenotypeTable' is a lazy view over the tables it was
# given. Copying it materializes the view, both because the view itself
# does not answer for the site scores the rest of TASSEL asks a table
# for, and because the copy keeps the sites of every table.
#
# The view holds its sites in the order the tables were handed over
# while reporting its position list sorted, and the copy inherits that
# pairing. Tables given out of genomic order, or holding interleaved
# sites, would therefore report one site's position against another
# site's calls. The sort that put the position list in order also
# reports where each sorted site came from, and filtering the copy
# through that redirect puts the two back in step.
#
# @param jGts
# A 'list' of Java 'GenotypeTable' objects.
# @param how
# The join method the caller was asked for: '"intersectJoin"' or
# '"unionJoin"'.
#
# @return
# A Java 'GenotypeTable' object.
.joinGenotypeTables <- function(jGts, how) {
    if (length(jGts) == 1L) {
        return(jGts[[1L]])
    }

    if (!all(vapply(jGts, function(jGt) jGt$hasGenotype(), logical(1)))) {
        rlang::abort(c(
            "Only genotype tables holding allele calls can be joined",
            "x" = paste0(
                "At least one table holds site scores - reference ",
                "probabilities or dosages - rather than allele calls."
            )
        ))
    }

    jArray <- rJava::.jarray(
        x              = jGts,
        contents.class = TASSEL_JVM$GENOTYPE_TABLE
    )

    # A Java 'Throwable' cannot be used as 'parent': rJava's '$' method
    # intercepts the fields rlang probes when formatting a chained error
    joinFailed <- function(cnd) {
        rlang::abort(c(
            "Could not join the genotype tables",
            "i" = "Do the tables share any taxa IDs?",
            "i" = paste0("TASSEL reported: ", conditionMessage(cnd))
        ))
    }

    jCombined <- tryCatch(
        rJava::J(TASSEL_JVM$COMBINE_GENOTYPE_TABLE)$getInstance(
            jArray,
            how == "unionJoin"
        ),
        error = joinFailed
    )

    # An intersect that keeps no taxa surfaces as an empty table rather
    # than an exception, so both outcomes report the likely cause
    if (jCombined$numberOfTaxa() == 0L) {
        rlang::abort(c(
            "Could not join the genotype tables",
            "i" = "Do the tables share any taxa IDs?"
        ))
    }

    jGt <- rJava::J(TASSEL_JVM$GENOTYPE_TABLE_BUILDER)$
        getGenotypeCopyInstance(jCombined)

    plBuilder <- rJava::.jnew(TASSEL_JVM$POSITION_LIST_BUILDER)
    for (jGtIn in jGts) {
        plBuilder$addAll(jGtIn$positions())
    }

    sorted   <- plBuilder$buildWithSiteRedirect()
    redirect <- sorted$getY()

    if (identical(redirect, seq_along(redirect) - 1L)) {
        return(jGt)
    }

    rJava::J(TASSEL_JVM$FILTER_GENOTYPE_TABLE)$getInstance(
        jGt,
        rJava::.jcast(sorted$getX(), TASSEL_JVM$POSITION_LIST),
        redirect
    )
}


## ----
# Join a collection of Java Phenotype and GenotypeTable objects
#
# @param coll
# The list returned by '.collectPhenotypes()'.
# @param how
# The 'PhenotypeBuilder' join method to call: '"intersectJoin"',
# '"unionJoin"', or '"concatenate"'.
#
# @return
# A 'TasselPhenotype', a 'TasselGenotype' when only genotype data was
# collected, a 'TasselGenomicDataset' when both were, or a
# 'TasselGenotypePhenotype' if every data-bearing input was one.
.joinPhenotypes <- function(coll, how) {
    jGt <- if (length(coll$jGts) == 0L) {
        NULL
    } else {
        .joinGenotypeTables(coll$jGts, how)
    }

    if (coll$jPhenotypes$isEmpty()) {
        if (coll$legacy) {
            return(.tasselObjectConstructor(jGt))
        }

        return(createTasselGenotype(jGt))
    }

    jPh <- if (coll$jPhenotypes$size() == 1L) {
        # TASSEL will not join a list of one, and there is nothing to join
        rJava::.jcast(coll$jPhenotypes$get(0L), TASSEL_JVM$PHENOTYPE)
    } else {
        builder <- rJava::.jnew(TASSEL_JVM$PHENO_BUILDER)$
            fromPhenotypeList(coll$jPhenotypes)

        builder <- switch(
            how,
            "intersectJoin" = builder$intersectJoin(),
            "unionJoin"     = builder$unionJoin(),
            "concatenate"   = builder$concatenate()
        )

        builder$build()$get(0L)
    }

    if (is.null(jGt)) {
        if (coll$legacy) {
            return(.tasselObjectConstructor(jPh))
        }

        return(createTasselPhenotype(jPh))
    }

    jGp <- joinGenotypePhenotype(
        jGt,
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
#' @title Intersect join phenotype and genotype tables
#'
#' @description Intersect join any number of phenotype objects based on the
#'    \code{Taxa} column. Objects carrying only genotype data are joined into
#'    a single genotype table holding the sites of each, and the taxa they
#'    all hold. If phenotype data was given as well, the joined phenotype is
#'    attached to that table and a
#'    \code{\linkS4class{TasselGenomicDataset}} is returned.
#'
#' @details Joining genotype tables is the way back from data split by
#'    chromosome or by collection of sites. Sites are returned in genomic
#'    order no matter which order the tables were given in, and the tables
#'    are expected to hold different sites - use
#'    \code{\link{mergeGenotypeTables}} to merge the calls of tables that
#'    describe the same sites. Depth and the other per-site scores are not
#'    carried into the joined table.
#'
#' @param ... Any number of rTASSEL objects containing phenotype or genotype
#'    data. Accepted classes are \code{\linkS4class{TasselPhenotype}},
#'    \code{\linkS4class{TasselGenotype}},
#'    \code{\linkS4class{TasselGenomicDataset}},
#'    \code{\linkS4class{PCAResults}}, \code{\linkS4class{MDSResults}},
#'    \code{\linkS4class{AssociationResultsBLUE}}, and the deprecated
#'    \code{TasselGenotypePhenotype}. Lists of objects are flattened, so
#'    earlier \code{intersectJoin(c(ph1, ph2))} style calls keep working.
#'
#' @return A \code{\linkS4class{TasselPhenotype}} object, a
#'    \code{\linkS4class{TasselGenotype}} object if only genotype data was
#'    given, or a \code{\linkS4class{TasselGenomicDataset}} if both were.
#'    Returns a \code{TasselGenotypePhenotype} if every data-bearing input
#'    was one.
#'
#' @seealso \code{\link{mergeGenotypeTables}}
#'
#' @examples
#' \dontrun{
#' # Merge several phenotype tables of covariates and traits
#' intersectJoin(ph1Cov, ph2Traits, ph3MoreTraits)
#'
#' # Put genotype data split by chromosome back together
#' intersectJoin(gtChr1, gtChr2, gtChr3)
#'
#' # Attach a genotype table at the same time
#' intersectJoin(gt, ph1Cov, ph2Traits, ph3MoreTraits)
#'
#' # Carry BLUE estimates forward as the traits of a new data set
#' blues <- assocModelFitter(ph2Traits, . ~ .)
#' intersectJoin(gt, blues, ph1Cov)
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
#' @title Union join phenotype and genotype tables
#'
#' @description Union join any number of phenotype objects based on the
#'    \code{Taxa} column. Objects carrying only genotype data are joined into
#'    a single genotype table holding the sites of each, and every taxon any
#'    of them holds. If phenotype data was given as well, the joined
#'    phenotype is attached to that table and a
#'    \code{\linkS4class{TasselGenomicDataset}} is returned.
#'
#' @details Joining genotype tables is the way back from data split by
#'    chromosome or by collection of sites. Sites are returned in genomic
#'    order no matter which order the tables were given in, and calls a
#'    table never made - those of a taxon another table alone holds - are
#'    returned as missing. The tables are expected to hold different sites -
#'    use \code{\link{mergeGenotypeTables}} to merge the calls of tables
#'    that describe the same sites. Depth and the other per-site scores are
#'    not carried into the joined table.
#'
#' @param ... Any number of rTASSEL objects containing phenotype or genotype
#'    data. Accepted classes are \code{\linkS4class{TasselPhenotype}},
#'    \code{\linkS4class{TasselGenotype}},
#'    \code{\linkS4class{TasselGenomicDataset}},
#'    \code{\linkS4class{PCAResults}}, \code{\linkS4class{MDSResults}},
#'    \code{\linkS4class{AssociationResultsBLUE}}, and the deprecated
#'    \code{TasselGenotypePhenotype}. Lists of objects are flattened, so
#'    earlier \code{unionJoin(c(ph1, ph2))} style calls keep working.
#'
#' @return A \code{\linkS4class{TasselPhenotype}} object, a
#'    \code{\linkS4class{TasselGenotype}} object if only genotype data was
#'    given, or a \code{\linkS4class{TasselGenomicDataset}} if both were.
#'    Returns a \code{TasselGenotypePhenotype} if every data-bearing input
#'    was one.
#'
#' @seealso \code{\link{mergeGenotypeTables}}
#'
#' @examples
#' \dontrun{
#' # Merge several phenotype tables of covariates and traits
#' unionJoin(ph1Cov, ph2Traits, ph3MoreTraits)
#'
#' # Put genotype data split by chromosome back together
#' unionJoin(gtChr1, gtChr2, gtChr3)
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
#'    \code{\linkS4class{PCAResults}}, \code{\linkS4class{MDSResults}},
#'    \code{\linkS4class{AssociationResultsBLUE}}, and the deprecated
#'    \code{TasselGenotypePhenotype}. Lists of objects are
#'    flattened, so earlier \code{concatenate(c(ph1, ph2))} style calls keep
#'    working.
#'    Unlike the joins, this function binds phenotype rows together and so
#'    does not accept genotype-only input. Genotype tables are combined by
#'    \code{\link{intersectJoin}} and \code{\link{unionJoin}}, which bring
#'    together the sites of each, or by
#'    \code{\link{mergeGenotypeTables}}, which merges the calls of tables
#'    describing the same sites.
#'
#' @return A \code{\linkS4class{TasselPhenotype}} object, or a
#'    \code{TasselGenotypePhenotype} if every data-bearing input was one.
#'
#' @seealso \code{\link{intersectJoin}}, \code{\link{unionJoin}}
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
#' @details
#' Tables are merged site by site, so a call two tables both make at the
#' same site is resolved into one. Use \code{\link{intersectJoin}} or
#' \code{\link{unionJoin}} instead to bring together tables that hold
#' different sites, such as data split by chromosome.
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
#' @seealso \code{\link{intersectJoin}}, \code{\link{unionJoin}}
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
