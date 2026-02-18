## ----
#' @title Wrapper function of TasselGenotypePhenotype class for genotype
#'    data
#'
#' @description This function is a wrapper for the
#'    \code{TasselGenotypePhenotype} class. It is used for storing genotype
#'    information into a class object.
#'
#' @name readGenotypeTableFromPath
#' @rdname readGenotypeTableFromPath
#'
#' @param path A genotype data path (e.g. \code{*.VCF, *.hmp}, etc.).
#' @param keepDepth Should depth be kept? Defaults to \code{FALSE}.
#' @param sortPositions Should positions be sorted? Defaults to \code{FALSE}.
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @importFrom rJava J
#' @importFrom rJava %instanceof%
#' @export
readGenotypeTableFromPath <- function(path, keepDepth = FALSE, sortPositions = FALSE) {
    warnMsg <- paste0(
        "The function 'readGenotypeTableFromPath()' will be deprecated soon.\n",
        "This will be replaced by '", cli::style_bold("readGenotype()"), "' in the next update."
    )
    message(warnMsg)

    if (!file.exists(path)) {
        stop("Cannot open file ", path, ": No such file or directory")
    }

    jrc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    return(
        .tasselObjectConstructor(
            jrc$read(path, keepDepth, sortPositions)
        )
    )
}


## ----
#' @title Create Summarized Experiment from a TASSEL Genotype Table
#'
#' @description This function will generate an object of
#'    \code{SummarizedExperiment} class for marker data derived from a
#'    \code{TasselGenotypePhenotype} class object.
#'
#' @name getSumExpFromGenotypeTable
#' @rdname getSumExpFromGenotypeTable
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#' @param coerceDosageToInt Should dosage array be returned as \code{integer}
#'    values? If \code{FALSE}, dosage array will be returned as type
#'    \code{raw} byte values. Returning \code{raw} byte values. Will greatly
#'    save on memory. Defaults to \code{TRUE}.
#' @param verbose Should messages be displayed to console? Defaults to
#'    \code{FALSE}.
#'
#' @return Returns a \code{SummarizedExperiment} of TASSEL genotype data.
#'
#' @importFrom rJava .jcall
#' @importFrom rJava is.jnull
#' @importFrom rJava .jevalArray
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
getSumExpFromGenotypeTable <- function(tasObj,
                                       coerceDosageToInt = TRUE,
                                       verbose = FALSE) {
    if (!inherits(tasObj, "TasselGenotypePhenotype")) {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGT <- getGenotypeTable(tasObj)

    if (rJava::is.jnull(jGT)) {
        stop("TASSEL genotype object not found")
    }

    # Create SumExp components (DF and ranges)
    sampleDF <- sampleDataFrame(tasObj)
    genomicRangesDF <- genomicRanges(jGT)

    # Create and return byte array from TASSEL
    if (verbose) message("Generating byte array...")
    jc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    genoCallByteArray <- jc$genotypeTableToDosageByteArray(jGT)
    if (verbose) message("Returning Java byte array to R...")
    dosMat <- lapply(genoCallByteArray, rJava::.jevalArray)
    if (coerceDosageToInt) {
        if (verbose) message("Coercing to integer...")
        dosMat <- lapply(dosMat, as.integer)

        # Replace 128 values (conversion artifact?) with NAs...
        dosMat <- lapply(dosMat, function(i) replace(i, i == 128, NA))
    }
    if (verbose) message("Transforming to SummarizedExperiment...")
    dosMat <- simplify2array(dosMat)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = dosMat,
        rowRanges = genomicRangesDF,
        colData = sampleDF
    )
    if (verbose) message("Finished.")
    warnMsg <- paste0("The function 'getSumExpFromGenotypeTable()' will be deprecated soon.")
    message(warnMsg, call. = FALSE)
    return(se)
}


## ----
## Get a GenotypeTable - not exported (house keeping)
getGenotypeTable <- function(jtsObject) {
    if(is(jtsObject, "TasselGenotypePhenotype")) {
        return(jtsObject@jGenotypeTable)
    }
    if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
    if(jtsObject %instanceof% "net.maizegenetics.dna.snp.GenotypeTable") {
        return(jtsObject)
    } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
        return(jtsObject$genotypeTable())
    } else {
        return(rJava::.jnull())
    }
}


## ----
## Return min/max physical positions from genotype tables (house keeping)
#' @importFrom rJava .jevalArray
#' @importFrom rJava is.jnull
getMinMaxPhysPositions <- function(tasObj) {
    if (!inherits(tasObj, "TasselGenotypePhenotype")) {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    javaGT <- getGenotypeTable(tasObj)

    positions <- javaGT$positions()
    chroms <- rJava::.jevalArray(javaGT$chromosomes())
    posLS  <- lapply(chroms, javaGT$firstLastSiteOfChromosome)

    physPos <- lapply(posLS, function(x) {
        firstPos <- positions$get(x[1])
        lastPos  <- positions$get(x[2])
        return(c(firstPos$getPosition(), lastPos$getPosition()))
    })
    names(physPos) <- sapply(chroms, function(x) x$getName())
    return(physPos)
}


## ----
## Return min/max physical positions from genotype tables (house keeping)
#' @importFrom rJava .jevalArray
#' @importFrom rJava is.jnull
getMinMaxVarSites <- function(tasObj) {
    if (!inherits(tasObj, "TasselGenotypePhenotype")) {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    javaGT <- getGenotypeTable(tasObj)

    positions <- javaGT$positions()
    chroms <- rJava::.jevalArray(javaGT$chromosomes())
    posLS  <- lapply(chroms, javaGT$firstLastSiteOfChromosome)

    names(posLS) <- sapply(chroms, function(x) x$getName())
    return(posLS)
}


##----
#' @title Read genotype data from GIGWA using QBMS
#'
#' @description Reads and stores genotype information from a
#'    \code{QBMS}-formatted data frame from a GIGWA server
#'
#' @param gigwa A \code{QBMS}-formatted GIGWA data frame object
#'
#' @importFrom rJava .jarray
#' @importFrom rJava J
#'
#' @export
readGenotypeTableFromGigwa <- function(gigwa) {
    warnMsg <- paste0(
        "The function 'readGenotypeTableFromGigwa()' will be deprecated soon.\n",
        "This will be replaced by '", cli::style_bold("readGenotype()"), "' in the next update."
    )
    message(warnMsg)

    plugin <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")

    matrixSub <- as.matrix(gigwa[, 5:ncol(gigwa)])
    mode(matrixSub) <- "integer"

    myGt <- plugin$createGenotypeFromRDataFrameElements(
        rJava::.jarray(colnames(gigwa[, 5:ncol(gigwa)])),
        rJava::.jarray(gigwa$chrom),
        rJava::.jarray(as.integer(gigwa$pos)),
        rJava::.jarray(gigwa$`rs#`),
        rJava::.jarray(gigwa$alleles),
        rJava::.jarray(matrixSub, dispatch = TRUE)
    )

    return(.tasselObjectConstructor(myGt))
}


## ----
#' @title Coerce genotype table to R matrix
#'
#' @description Converts a \code{TasselGenotypePhenotype} class into an R
#'    matrix if it contains genotype data.
#'
#' @param x A \code{TasselGenotypePhenotype} object
#' @param ... Additional arguments to be passed to or from methods.
#'
#' @importFrom rJava .jevalArray
#'
#' @export
as.matrix.TasselGenotypePhenotype <- function(x, ...) {
    plugin <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")

    if (!inherits(x, "TasselGenotypePhenotype")) {
        stop("`x` must be of class `TasselGenotypePhenotype`")
    }

    if (rJava::is.jnull(x@jGenotypeTable)) {
        stop("`x` must contain genotype data")
    }

    jg <- x@jGenotypeTable
    m <- rJava::.jevalArray(plugin$genotypeTableToDosageByteArray(jg), simplify = TRUE)
    mode(m) <- "integer"

    siteNames <- positionList(x)

    m[m == 128] <- NA
    colnames(m) <- siteNames$Name
    rownames(m) <- getTaxaIDs(x)

    return(m)
}


## ----
#' @title Get site summary of genotype table
#'
#' @description Returns positional data from a \code{TasselGenotypePhenotype}
#'    object
#'
#' @param tasObj A \code{TasselGenotypePhenotype} object
#'
#' @export
siteSummary <- function(tasObj) {
    if (!inherits(tasObj, "TasselGenotypePhenotype")) {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    if (rJava::is.jnull(tasObj@jGenotypeTable)) {
        stop("`tasObj` must contain genotype data")
    }

    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.data.GenotypeSummaryPlugin"),
        rJava::.jnull(),
        FALSE
    )

    plugin$setParameter("overview", tolower(as.character(FALSE)))
    plugin$setParameter("siteSummary", tolower(as.character(TRUE)))
    plugin$setParameter("taxaSummary", tolower(as.character(FALSE)))

    dataSet <- rJava::J("net.maizegenetics.plugindef.DataSet")
    summaryResults <- plugin$processData(dataSet$getDataSet(tasObj@jGenotypeTable))

    return(
        tableReportToDF(
            summaryResults$getData(0L)$getData()
        )
    )
}


## ----
#' @title Get taxa summary of genotype table
#'
#' @description Returns taxa data from a \code{TasselGenotypePhenotype}
#'    object
#'
#' @param tasObj A \code{TasselGenotypePhenotype} object
#'
#' @export
taxaSummary <- function(tasObj) {
    if (!inherits(tasObj, "TasselGenotypePhenotype")) {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    if (rJava::is.jnull(tasObj@jGenotypeTable)) {
        stop("`tasObj` must contain genotype data")
    }

    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.data.GenotypeSummaryPlugin"),
        rJava::.jnull(),
        FALSE
    )

    plugin$setParameter("overview", tolower(as.character(FALSE)))
    plugin$setParameter("siteSummary", tolower(as.character(FALSE)))
    plugin$setParameter("taxaSummary", tolower(as.character(TRUE)))

    dataSet <- rJava::J("net.maizegenetics.plugindef.DataSet")
    summaryResults <- plugin$processData(dataSet$getDataSet(tasObj@jGenotypeTable))

    return(
        tableReportToDF(
            summaryResults$getData(0L)$getData()
        )
    )
}


