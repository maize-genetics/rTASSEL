## ----
#' @title Wrapper function of TasselGenotypePhenotype class for
#'    GenotypePhenotype combined data
#'
#' @description
#' \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
#'
#' Creates a Java GenotypePhenotype object which is used for
#'    the deprecated \code{TasselGenotypePhenotype} object construction. The
#'    Java GenotypePhenotype object is created via an \code{intersect} method
#'    from TASSEL. Use \code{\link{readGenomicDataset}()} instead, which
#'    returns a \code{\linkS4class{TasselGenomicDataset}} object.
#'
#' @name readGenotypePhenotype
#' @rdname readGenotypePhenotype
#'
#' @param genoPathOrObj a path to a genotype file (e.g. VCF, hmp, etc.) or
#'    TASSEL Genotype Obj
#' @param phenoPathDFOrObj a path, a data frame of phenotypic data, or TASSEL
#'    Phenotype Obj
#' @param ... Additional parameters to be sent to the function. Currently,
#'    if an R data frame object is passed, additional parameters will also
#'    need to be entered for this process. These parameters are derived from
#'    the \code{readPhenotypeFromDataFrame} function. Mainly, \code{taxaID}
#'    is required. If you would like to specify depth retention and position
#'    sorting in Genotype Tables from a path, indicate them here. See
#'    \code{readGenotypeTableFromPath()} for more detail.
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @seealso \code{\link{readGenomicDataset}}
#'
#' @importFrom rJava J
#' @importFrom rJava is.jnull
#' @importFrom rJava .jinstanceof
#' @export
readGenotypePhenotype <- function(genoPathOrObj, phenoPathDFOrObj, ...) {
    lifecycle::deprecate_warn(
        "0.14.0", "readGenotypePhenotype()", "readGenomicDataset()"
    )

    # The legacy readers this delegates to are deprecated in their own right;
    # one notice per user-facing call is enough
    oldOpts <- options(lifecycle_verbosity = "quiet")
    on.exit(options(oldOpts), add = TRUE)

    genoObj <- getGenotypeTable(genoPathOrObj)
    if(rJava::is.jnull(genoObj)) {
        genoObj <- getGenotypeTable(
            readGenotypeTableFromPath(genoPathOrObj, keepDepth = FALSE, sortPositions = FALSE)
        )
    }
    phenoObj <- getPhenotypeTable(phenoPathDFOrObj)
    if(rJava::is.jnull(phenoObj) & is.data.frame(phenoPathDFOrObj)) {
        phenoObj <- getPhenotypeTable(
            readPhenotypeFromDataFrame(phenoPathDFOrObj, ...)
        )
    } else if(rJava::is.jnull(phenoObj) & !is.data.frame(phenoPathDFOrObj)) {
        phenoObj <- rJava::new(
            rJava::J("net/maizegenetics/phenotype/PhenotypeBuilder")
        )$fromFile(phenoPathDFOrObj)$build()$get(0L)
    }

    t <- .tasselObjectConstructor(
        combineTasselGenotypePhenotype(
            genotypeTable = genoObj,
            phenotype = phenoObj
        )
    )
    return(t)
}


## Combine TASSEL GenotypeTable and Phenotype to GenotypePhenotype class - not exported (house keeping)
combineTasselGenotypePhenotype <- function(genotypeTable, phenotype) {
    genotypeTable <- getGenotypeTable(genotypeTable)
    phenotype <- getPhenotypeTable(phenotype)
    new(J("net.maizegenetics.phenotype.GenotypePhenotypeBuilder"))$
        genotype(genotypeTable)$phenotype(phenotype)$intersect()$build()
}


## Get GenotypePhenotype - not exported - not exported (house keeping)
getGenotypePhenotype <- function(jtsObject) {
    if(is(jtsObject, "TasselGenotypePhenotype")) {
        jtsObject <- jtsObject@jTasselObj
    }
    jtsObject <- .unwrapTasselObject(jtsObject)
    if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
    if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
        return(jtsObject)
    } else {
        return(rJava::.jnull())
    }
}
