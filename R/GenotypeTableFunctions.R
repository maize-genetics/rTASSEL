#---------------------------------------------------------------------
# Script Name:   GenotypeTableFunctions.R
# Description:   Support working with TASSEL GenotypeTables
# Author:        Brandon Monier & Ed buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2022-04-28 at 18:02:36
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions
#    necessary for reading in genotype datasets into R and
#    extracting data from TASSEL genotype objects.
#--------------------------------------------------------------------

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



#' @title Read \code{GenotypeTable} objects from BrAPI
#'
#' @description Reads \code{GenotypeTable} objects from BrAPI
#'
#' @name readGenotypeTableFromBrapi
#' @rdname readGenotypeTableFromBrapi
#'
#' @param brapiObj A \code{BrapiCon} object.
#' @param variantSetDbId What set of variants do you want to read from BrAPI?
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @export
readGenotypeTableFromBrapi <- function(
    brapiObj,
    variantSetDbId = NULL
) {
    if (class(brapiObj) != "BrapiCon") {
        stop("`brapiObj` must be of class `BrapiCon`")
    }

    if (is.null(variantSetDbId)) {
        stop("You must choose a BrAPI database ID")
    }

    ## Get JSON components ----
    baseUrl <- rPHG::brapiURL(brapiObj)
    variantMatUrl <- paste0(baseUrl, "/variantmatrix?variantSetDbId=", variantSetDbId)
    vMJson <- rPHG:::parseJSON(variantMatUrl)

    ## Get BrAPI components ----
    ## TODO - how to deal with phasing?
    snpRMat <- gsub("/|\\|", "", t(vMJson$result$data))
    taxa    <- vMJson$result$callSetDbIds
    posDf   <- vMJson$result$variants


    ## Internal parameters ----
    ## TODO - how to get actual alleles for each site?
    ## Defaults to A, C, or N
    alleleStates <- list(
        ref  = "A",
        alt  = "C",
        null = "N"
    )

    ## Instantiate TASSEL API classes ----
    nucConst         <- rJava::J("net/maizegenetics/dna/snp/NucleotideAlignmentConstants")
    gtUtils          <- rJava::J("net/maizegenetics/dna/snp/GenotypeTableUtils")
    genoCallBuilder  <- rJava::J("net/maizegenetics/dna/snp/genotypecall/GenotypeCallTableBuilder")
    genoTableBuilder <- rJava::J("net/maizegenetics/dna/snp/GenotypeTableBuilder")
    posListBuilder   <- rJava::.jnew("net/maizegenetics/dna/map/PositionListBuilder")
    taxaListBuilder  <- rJava::.jnew("net/maizegenetics/taxa/TaxaListBuilder")


    ## GenotypeCallTable ----
    genotypeCalls <- genoCallBuilder$getUnphasedNucleotideGenotypeBuilder(
        as.integer(nrow(snpRMat)),
        as.integer(ncol(snpRMat))
    )
    for (i in seq_len(nrow(snpRMat))) {
        for (j in seq_len(ncol(snpRMat))) {

            # TODO - how to deal with more than one alt state?
            curAllele <- switch (
                EXPR = snpRMat[i, j],
                "00" = c(alleleStates$ref, alleleStates$ref),
                "10" = c(alleleStates$alt, alleleStates$ref),
                "01" = c(alleleStates$ref, alleleStates$alt),
                "11" = c(alleleStates$alt, alleleStates$alt),
                `NA` = c(alleleStates$null, alleleStates$null)
            )

            byteVal <- gtUtils$getDiploidValue(
                rJava::.jbyte(nucConst$getNucleotideDiploidByte(curAllele[1])),
                rJava::.jbyte(nucConst$getNucleotideDiploidByte(curAllele[2]))
            )
            genotypeCalls$setBase(
                as.integer(i - 1),
                as.integer(j - 1),
                rJava::.jbyte(byteVal)
            )
        }
    }


    ## TaxaList ----
    for (i in taxa) {
        taxaListBuilder$add(rJava::.jnew("net/maizegenetics/taxa/Taxon", i))
    }


    ## PositionList ----
    for (i in seq_len(nrow(posDf))) {
        posListBuilder$add(
            rJava::J(
                "net/maizegenetics/dna/map/Position", "of",
                posDf[i, ]$contig,
                as.integer(posDf[i, ]$start)
            )
        )
    }


    ## Build CoreGenotypeTable ----
    myGt <- genoTableBuilder$getInstance(
        genotypeCalls$build(),
        posListBuilder$build(),
        taxaListBuilder$build()
    )

    ## Return `TasselGenotypePhenotype` object (rTASSEL) ----
    return(.tasselObjectConstructor(myGt))
}


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
#' @param useRef Use reference alleles or major alleles at sites. If
#'    \code{FALSE}, major alleles will be used. Defaults to \code{FALSE}.
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
                                       useRef = FALSE,
                                       coerceDosageToInt = TRUE,
                                       verbose = FALSE) {

    if (class(tasObj) != "TasselGenotypePhenotype") {
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
    genoCallByteArray <- jc$genotypeTableToDosageByteArray(
        jGT,
        useRef
    )
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
    return(se)
}


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


## Return min/max physical positions from genotype tables (house keeping)
#' @importFrom rJava .jevalArray
#' @importFrom rJava is.jnull
getMinMaxPhysPositions <- function(tasObj) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
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


## Return min/max physical positions from genotype tables (house keeping)
#' @importFrom rJava .jevalArray
#' @importFrom rJava is.jnull
getMinMaxVarSites <- function(tasObj) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
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


