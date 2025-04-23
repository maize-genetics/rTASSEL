## ----
#' @title Imputation methods in Numerical Transformations
#'
#' @description This method takes an input \code{GenotypeTable} object with
#'    missing values and imputes the missing values using one of the chosen
#'    methods.
#'
#' @param tasObj an rTASSEL \code{TasselGenotypePhenotype} object.
#' @param byMean Will imputation be performed by computing the mean of the
#'    respective column? Defaults to \code{TRUE}.
#' @param nearestNeighbors Number of nearest neighbors to be evaluated.
#'    Defaults to \code{5}.
#' @param distance Distance type. Options are \code{Euclidean},
#'    \code{Manhattan}, or \code{Cosine}. Defaults to \code{Euclidean}.
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
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

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

    jPhenoTable <- getPhenotypeTable(tasObj)
    if (rJava::is.jnull(jPhenoTable)) {
        return(.tasselObjectConstructor(res$getData(0L)$getData()))
    } else {
        return(
            .tasselObjectConstructor(
                combineTasselGenotypePhenotype(
                    genotypeTable = res$getData(0L)$getData(),
                    phenotype = jPhenoTable
                )
            )
        )
    }
}



#' @title LD KNNi imputation
#'
#' @description This imputation algorithm uses LD to identify good predictors
#'    for each SNP, and then uses the high LD SNPs to identify K-Nearest
#'    Neighbors. The genotype is called with a weighted mode of the KNNs.
#'
#' @param tasObj an rTASSEL \code{TasselGenotypePhenotype} object.
#' @param highLDSSites Number of sites in high LD to use in imputation.
#'    Acceptable values are between \code{2} and \code{2000}. Defaults to
#'    \code{30}.
#' @param knnTaxa Number of neighbors to use in imputation. Acceptable values
#'    are between \code{2} and \code{200}. Defaults to \code{10}.
#' @param maxDistance Maximum physical distance between sites to search for LD
#'    (-1 for no distance cutoff - unlinked chromosomes will be tested).
#'    Defaults to \code{10e6}.
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
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

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

    jPhenoTable <- getPhenotypeTable(tasObj)
    if (rJava::is.jnull(jPhenoTable)) {
        return(.tasselObjectConstructor(res$getData(0L)$getData()))
    } else {
        return(
            .tasselObjectConstructor(
                combineTasselGenotypePhenotype(
                    genotypeTable = res$getData(0L)$getData(),
                    phenotype = jPhenoTable
                )
            )
        )
    }
}


