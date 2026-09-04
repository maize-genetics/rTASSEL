## ----
#' @title Imputation methods in Numerical Transformations
#'
#' @description This method takes an input \code{GenotypeTable} object with
#'    missing values and imputes the missing values using one of the chosen
#'    methods.
#'
#' @param tasObj An object of class \code{\linkS4class{TasselGenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the deprecated
#'    \code{TasselGenotypePhenotype} class are still accepted.
#' @param byMean Will imputation be performed by computing the mean of the
#'    respective column? Defaults to \code{TRUE}.
#' @param nearestNeighbors Number of nearest neighbors to be evaluated.
#'    Defaults to \code{5}.
#' @param distance Distance type. Options are \code{Euclidean},
#'    \code{Manhattan}, or \code{Cosine}. Defaults to \code{Euclidean}.
#'
#' @return An object of the same class as \code{tasObj} holding the imputed
#'    genotype data.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#'
#' @export
imputeNumeric <- function(
    tasObj,
    byMean = TRUE,
    nearestNeighbors = 5,
    distance = c("Euclidean", "Manhattan", "Cosine")
) {
    tasIn      <- .resolveTasselInput(tasObj, "genotype", "imputeNumeric")
    jGenoTable <- tasIn$jGt

    distance <- match.arg(distance)

    # Create imputation plugin
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.numericaltransform.ImputationPlugin"),
        rJava::.jnull(),
        FALSE
    )

    # Set parameters
    plugin$setParameter("ByMean", tolower(as.character(byMean)))
    plugin$setParameter("nearestNeighbors", as.character(nearestNeighbors))
    plugin$setParameter("distance", distance)

    dataSet <- rJava::J("net.maizegenetics.plugindef.DataSet")
    res <- plugin$processData(dataSet$getDataSet(jGenoTable))

    .wrapImputed(res$getData(0L)$getData(), tasIn)
}


## ----
# Rebuild an imputation result as the class that was handed in
#
# @description
# Imputation plugins hand back a bare `GenotypeTable`. Inputs that also
# carried phenotype data are re-joined against it so that the returned
# object holds the same components as the original.
#
# @param jGt
# The Java `GenotypeTable` produced by the imputation plugin.
# @param tasIn
# The list returned by `.resolveTasselInput()` for the original input.
#
# @return
# An object of the same class as `tasIn$original`.
.wrapImputed <- function(jGt, tasIn) {
    jRes <- if (rJava::is.jnull(tasIn$jPh)) {
        jGt
    } else {
        combineTasselGenotypePhenotype(
            genotypeTable = jGt,
            phenotype     = tasIn$jPh
        )
    }

    .wrapLikeInput(jRes, tasIn$original)
}



#' @title LD KNNi imputation
#'
#' @description This imputation algorithm uses LD to identify good predictors
#'    for each SNP, and then uses the high LD SNPs to identify K-Nearest
#'    Neighbors. The genotype is called with a weighted mode of the KNNs.
#'
#' @param tasObj An object of class \code{\linkS4class{TasselGenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the deprecated
#'    \code{TasselGenotypePhenotype} class are still accepted.
#' @param highLDSSites Number of sites in high LD to use in imputation.
#'    Acceptable values are between \code{2} and \code{2000}. Defaults to
#'    \code{30}.
#' @param knnTaxa Number of neighbors to use in imputation. Acceptable values
#'    are between \code{2} and \code{200}. Defaults to \code{10}.
#' @param maxDistance Maximum physical distance between sites to search for LD
#'    (-1 for no distance cutoff - unlinked chromosomes will be tested).
#'    Defaults to \code{10e6}.
#'
#' @return An object of the same class as \code{tasObj} holding the imputed
#'    genotype data.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#'
#' @export
imputeLDKNNi <- function(
    tasObj,
    highLDSSites = 30,
    knnTaxa = 10,
    maxDistance = 10e6
) {
    tasIn      <- .resolveTasselInput(tasObj, "genotype", "imputeLDKNNi")
    jGenoTable <- tasIn$jGt

    if (highLDSSites < 2 || highLDSSites > 2000) {
        stop("`highLDSSites` exceeds acceptable parameters (2 - 2000).")
    }
    if (knnTaxa < 2 || knnTaxa > 200) {
        stop("`knnTaxa` exceeds acceptable parameters (2 - 200).")
    }

    # Create imputation plugin
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.imputation.LDKNNiImputationPlugin"),
        rJava::.jnull(),
        FALSE
    )

    # Set parameters
    plugin$setParameter("highLDSSites", as.character(highLDSSites))
    plugin$setParameter("knnTaxa", as.character(knnTaxa))
    plugin$setParameter("maxLDDistance", as.character(as.integer(maxDistance)))

    dataSet <- rJava::J("net.maizegenetics.plugindef.DataSet")
    res <- plugin$processData(dataSet$getDataSet(jGenoTable))

    .wrapImputed(res$getData(0L)$getData(), tasIn)
}


