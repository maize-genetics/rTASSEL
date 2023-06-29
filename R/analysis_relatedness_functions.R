#--------------------------------------------------------------------
# Script Name:   AnalysisRelatednessFunctions.R
# Description:   Functions for TASSEL relatedness analyses
# Author:        Brandon Monier
# Created:       2019-04-04 at 21:31:09
# Last Modified: 2022-03-01 at 17:50:19
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions necessary
#    for relatedness-based methods in TASSEL.
#--------------------------------------------------------------------

#' @title Create a TASSEL kinship matrix
#'
#' @description This function will calculate a kinship matrix using
#'    TASSEL's \code{Kinship} plugin and respective parameters.
#'
#' @name kinshipMatrix
#' @rdname kinshipMatrix
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#' @param method A Kinship method. Defaults to \code{Centered_IBS}. Other
#'    options include \code{Normalized_IBS}, \code{Dominance_Centered_IBS},
#'    and \code{Dominance_Normalized_IBS}.
#' @param maxAlleles Maximum number of alleles. Can be within the range of
#'    \code{2} to \code{6}.
#' @param algorithmVariation Algorithm variation. If
#'    \code{Dominance_Centered_IBS} is selected, users can switch between
#'    \code{Observed_Allele_Freq} and \code{Proportion_Heterozygous}.
#'
#' @return Returns a `TasselDistanceMatrix` object.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#' @export
kinshipMatrix <- function(tasObj,
                          method = "Centered_IBS",
                          maxAlleles = 6,
                          algorithmVariation = "Observed_Allele_Freq") {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    # Create kinship plugin
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.distance.KinshipPlugin"),
        rJava::.jnull(),
        FALSE
    )
    plugin$setParameter("method", toString(method))
    plugin$setParameter("maxAlleles", toString(maxAlleles))
    plugin$setParameter("algorithmVariation", toString(algorithmVariation))
    distMatrix <- plugin$runPlugin(jGenoTable)

    jTl <- distMatrix$getTaxaList()

    tl <- sapply(1:distMatrix$numberOfTaxa(), function(i) {
        jTl$taxaName(as.integer(i - 1))
    })

    methods::new(
        Class = "TasselDistanceMatrix",
        taxa = tl,
        numTaxa = distMatrix$numberOfTaxa(),
        summaryMatrix = summaryDistance(distMatrix),
        jDistMatrix = distMatrix
    )
}


#' @title Create a TASSEL distance matrix
#'
#' @description This function will calculate a distance matrix using
#'    TASSEL's \code{DistanceMatrix} plugin.
#'
#' @name distanceMatrix
#' @rdname distanceMatrix
#'
#' @param tasObj an rTASSEL \code{TasselGenotypePhenotype} object
#'
#' @return Returns a `TasselDistanceMatrix` object.
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#' @export
distanceMatrix <- function(tasObj) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.distance.DistanceMatrixPlugin"),
        rJava::.jnull(),
        FALSE
    )
    distMatrix <- plugin$getDistanceMatrix(jGenoTable)

    jTl <- distMatrix$getTaxaList()

    tl <- sapply(1:distMatrix$numberOfTaxa(), function(i) {
        jTl$taxaName(as.integer(i - 1))
    })

    methods::new(
        Class = "TasselDistanceMatrix",
        taxa = tl,
        numTaxa = distMatrix$numberOfTaxa(),
        summaryMatrix = summaryDistance(distMatrix),
        jDistMatrix = distMatrix
    )
}


#' @title read TASSEL distance matrix object from file
#'
#' @description This function will read a TASSEL distance matrix from
#'    file and convert it into a \code{TasselDistanceMatrix} object.
#'
#' @name readTasselDistanceMatrix
#' @rdname readTasselDistanceMatrix
#'
#' @param file A file path of type \code{character}
#'
#' @return Returns a \code{TasselDistanceMatrix} object.
#'
#' @importFrom methods new
#' @importFrom rJava J
#'
#' @export
readTasselDistanceMatrix <- function(file) {

    if (!file.exists(file)) {
        stop("File does not exist.", call. = FALSE)
    }

    rJC <- rJava::J("net/maizegenetics/taxa/distance/ReadDistanceMatrix")
    distMatrix <- rJC$readDistanceMatrix(file)

    jTl <- distMatrix$getTaxaList()

    tl <- sapply(1:distMatrix$numberOfTaxa(), function(i) {
        jTl$taxaName(as.integer(i - 1))
    })

    methods::new(
        Class = "TasselDistanceMatrix",
        taxa = tl,
        numTaxa = distMatrix$numberOfTaxa(),
        summaryMatrix = summaryDistance(distMatrix),
        jDistMatrix = distMatrix
    )
}


#' @title Coerce matrix to TasselDistanceMatrix object
#'
#' @description Coerces an object of \code{matrix} class into an rTASSEL
#'    object of \code{TasselDistanceMatrix} class.
#'
#' @param m An object of \code{matrix} class of \eqn{m \times m} structure
#'    (e.g. a pairwise matrix). Additionally, row and column names must be
#'    the same.
#'
#' @return Returns a \code{TasselDistanceMatrix} object.
#'
#' @importFrom rJava .jarray
#' @importFrom rJava J
#' @importFrom methods new
#'
#' @export
asTasselDistanceMatrix <- function(m) {

    if (!inherits(m, "matrix")) {
        stop("'m' parameter must be a matrix object", call. = FALSE)
    }
    if (nrow(m) != ncol(m)) {
        stop("Matrix object must have equal rows and columns", call. = FALSE)
    }
    if (is.null(colnames(m)) || is.null(rownames(m))) {
        stop("Matrix object must have column and row names", call. = FALSE)
    }
    if (!all(colnames(m) == rownames(m))) {
        stop("Matrix object must have the same row and column name structure", call. = FALSE)
    }

    plugin <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")

    mJ <- rJava::.jarray(m, dispatch = TRUE)
    tJ <- rJava::.jarray(colnames(m), dispatch = TRUE)

    distMatrix <- plugin$asTasselDistanceMatrix(mJ, tJ)

    methods::new(
        Class = "TasselDistanceMatrix",
        taxa = colnames(m),
        numTaxa = distMatrix$numberOfTaxa(),
        summaryMatrix = summaryDistance(distMatrix),
        jDistMatrix = distMatrix
    )

}


#' @title Run PCA on Genotype Table
#'
#' @description This method performs principal components analysis and returns
#'    the requested number of PC axes (components).
#'
#' @param tasObj an rTASSEL \code{TasselGenotypePhenotype} object.
#' @param useCovariance If \code{TRUE}, analysis will do an eigenvalue
#'    decomposition of the covariance matrix. If \code{FALSE}, it will use a
#'    correlation matrix. NOTE: Using the covariance matrix is recommended for
#'    genotypes while the correlation matrix is often used for phenotypes.
#'    Defaults to \code{TRUE}.
#' @param limitBy This parameter determines the type of value that will be used
#'    to limit the number of principal components (axes) returned. The possible
#'    choices are \code{number_of_components}, \code{min_eigenvalue},
#'    and \code{total_variance}.
#' @param nComponents The analysis will return this many principal components
#'    up to the number of taxa.
#' @param minEigenval All principal components with an eigenvalue greater than
#'    or equal to this value will be returned. NOTE: works only if
#'    \code{min_eigenvalue} is set in the \code{limitBy} parameter.
#' @param totalVar The first principal components that together explain this
#'    proportion of the total variance will be returned. NOTE: works only if
#'    \code{total_variance} is set in the \code{limitBy} parameter.
#' @param reportEigenvalues Returns a list of eigenvalues sorted high to low.
#' @param reportEigenvectors Returns the eigenvectors calculated from a
#'    Singular Value Decomposition of the data. The resulting table can be
#'    quite large if the number of variants and taxa are big.
#'
#' @return A \code{DataFrame} object.
#'
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#'
#' @export
pca <- function(
    tasObj,
    useCovariance = TRUE,
    limitBy = c("number_of_components", "min_eigenvalue", "total_variance"),
    nComponents = 5,
    minEigenval = 0,
    totalVar = 0.5,
    reportEigenvalues = TRUE,
    reportEigenvectors = TRUE
) {

    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    limitBy <- match.arg(limitBy)

    # Create PCA plugin
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.data.PrincipalComponentsPlugin"),
        rJava::.jnull(),
        FALSE
    )

    # Set PCA parameters
    plugin$setParameter("covariance", tolower(as.character(useCovariance)))
    plugin$setParameter("limitBy", limitBy)
    plugin$setParameter("ncomponents", as.character(nComponents))
    plugin$setParameter("minEigenval", as.character(minEigenval))
    plugin$setParameter("reportEigenvalues", tolower(as.character(reportEigenvalues)))
    plugin$setParameter("reportEigenvectors", tolower(as.character(reportEigenvectors)))

    # Run PCA plugin
    dataSet <- rJava::J("net.maizegenetics.plugindef.DataSet")
    pcaRes <- plugin$processData(dataSet$getDataSet(jGenoTable))

    reportBody <- lapply(seq_len(pcaRes$getSize()), function(i) {
        b <- pcaRes$getData(as.integer(i - 1))
        return(tableReportToDF(b$getData()))
    })

    reportNames <- lapply(seq_len(pcaRes$getSize()), function(i) {
        pcaRes$getData(as.integer(i - 1))$getName()
    })

    names(reportBody) <- unlist(reportNames)

    if (reportEigenvalues) {
        colnames(reportBody$Eigenvalues_Datum) <- gsub(" ", "_", colnames(reportBody$Eigenvalues_Datum))
    }



    # Add `CorePhenotype` object to return
    reportBody$jPcaObj <- pcaRes$getDataWithName("PC_Datum")$get(0L)$getData()

    if (!reportEigenvalues && !reportEigenvectors) {
        return(reportBody[c("PC_Datum", "jPcaObj")])
    } else {
        return(reportBody)
    }
}


#' @title Run MDS on \code{TasselDistanceMatrix} objects
#'
#' @description Perform multidimensional scaling (MDS) on
#'    \code{TasselDistanceMatrix} objects.
#'
#' @param distMat A \code{TasselDistanceMatrix} object.
#' @param nAxes The number of axes or dimensions and associated eigenvalues to
#'    be returned by the analysis. Defaults to \code{5}.
#' @param removeNaN Remove \code{NaNs} from matrix before performing MDS.
#'    Defaults to \code{TRUE}.
#'
#' @return A \code{DataFrame} object.
#'
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#'
#' @export
mds <- function(
    distMat,
    nAxes = 5,
    removeNaN = TRUE
) {
    if (class(distMat) != "TasselDistanceMatrix") {
        stop("`distMat` must be of class `TasselDistanceMatrix`")
    }

    # Create MDS plugin
    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.distance.MultiDimensionalScalingPlugin"),
        rJava::.jnull(),
        FALSE
    )

    # Set parameters
    plugin$setParameter("axes", as.character(nAxes))
    plugin$setParameter("removeNaN", tolower(as.character(removeNaN)))

    # Run PCA plugin
    dataSet <- rJava::J("net.maizegenetics.plugindef.DataSet")
    mdsRes <- plugin$performFunction(dataSet$getDataSet(distMat@jDistMatrix))

    return(tableReportToDF(mdsRes$getData(0L)$getData()))
}


