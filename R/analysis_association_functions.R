## ----
#' @title R interface for TASSEL's association methods
#'
#' @description This function acts as a front-end for TASSEL's extensive
#'    association analysis methods. Using this function, users can run
#'    the following TASSEL association methods:
#'    \itemize{
#'      \item{best linear unbiased estimates (BLUEs)}
#'      \item{generalized linear model (GLM)}
#'      \item{mixed linear model}
#'      \item{Fast association (Shabalin 2012)}
#'    }
#'
#' @name assocModelFitter
#' @rdname assocModelFitter
#'
#' @param tasObj An object of class \code{\linkS4class{TasselPhenotype}} or
#'   \code{\linkS4class{TasselGenomicDataset}}. A genotype table is required
#'   for every analysis except BLUEs. Objects of the deprecated
#'   \code{TasselGenotypePhenotype} class are still accepted.
#' @param formula An R-based linear model formula. The general layout of this
#'   formula uses the following TASSEL data scheme:
#'   \code{<data> ~ <factor> and/or <covariate>}. If all traits in a Phenotype
#'   object should be ran, a simplified formula (\code{. ~ .}) can be used.
#'   This scheme can also be used for running all \code{<data>} or
#'   \code{<factor>} and/or \code{<covariate>} data as well. Single variables
#'   are separated witha \code{+} operator. See vignette for further
#'   clarification.
#' @param fitMarkers Should marker data be fitted? If \code{TRUE}, GLM
#'   analysis will be executed. If \code{FALSE}, BLUEs will be calculated.
#'   Defaults to \code{FALSE}.
#' @param kinship Should kinship data be accounted for in the model? If so,
#'   a TASSEL kinship matrix object of class \code{TasselDistanceMatrix} must
#'   be submitted. Defaults to \code{NULL}.
#' @param fastAssociation Should TASSEL's Fast Association plugin be used?
#'   Consider setting to \code{TRUE} if you have many phenotypes in your
#'   data set.
#' @param maxP Maximum p-value (0 - 1) to be reported. Currently works with
#'   fast association only. Defaults to a p-value of \code{0.001}
#'   will be used as a threshold. \strong{Note:} p-value parameter will
#'   not be used for BLUE analysis.
#' @param maxThreads Maximum threads to be used when running fast association.
#'   If \code{NULL}, all threads on machine will be used.
#' @param outputFile Output file prefix to be specified in case you want
#'   to write data directly to disk. Highly recommended for large datasets.
#'   If \code{NULL}, no data will be saved to disk. If a character
#' @param minClassSize The minimum acceptable genotype class size. Genotypes
#'   in a class with a smaller size will be set to missing. Defaults to 0.
#' @param biallelicOnly Only test sites that are bi-allelic. The alternative is
#'   to test sites with two or more alleles. Defaults to \code{FALSE}
#' @param appendAddDom If true, additive and dominance effect estimates will
#'   be added to the stats report for bi-allelic sites only. The effect will
#'   only be estimated when the data source is genotype (not a probability).
#'   The additive effect will always be non-negative. Defaults to \code{FALSE}.
#'
#' @return Returns an R list containing \code{DataFrame}-based data frames
#'
#' @importFrom rJava is.jnull
#' @importFrom rJava J
#' @importFrom rJava .jnew
#' @importFrom rJava .jnull
#' @importFrom rJava new
#' @importFrom rlang .data
#' @export
assocModelFitter <- function(
    tasObj,
    formula,
    fitMarkers = FALSE,
    kinship = NULL,
    fastAssociation = FALSE,
    maxP = 0.001,
    maxThreads = 1,
    minClassSize = 0,
    outputFile = NULL,
    biallelicOnly = FALSE,
    appendAddDom = FALSE
) {

    tasIn <- .resolveTasselInput(tasObj, "phenotype", "assocModelFitter")
    jGt   <- tasIn$jGt
    jPh   <- tasIn$jPh

    # Logic - check kinship object
    if (!is.null(kinship) && !inherits(kinship, "TasselDistanceMatrix")) {
        stop("TASSEL kinship object is not of TasselDistanceMatrix class", call. = FALSE)
    }
    if (!is.null(kinship) && inherits(kinship, "TasselDistanceMatrix")) {
        kinTaxa <- colnames(kinship)
        genoTaxa <- getTaxaIDs(tasObj)

        if (!any(kinTaxa %in% genoTaxa)) {
            stop("No taxa IDs in your kinship object match your genotype information.", call. = FALSE)
        } else {
            kinship <- kinship@jDistMatrix
        }
    }

    # Subset phenotype data
    rData        <- tableReportToDF(jPh)
    attrData     <- makeAttributeData(jPh, rData)
    traitsToKeep <- unlist(parseFormula(formula, attrData))

    # Logic - Handle association analyses
    jRC <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    jTasFilt <- tasPhenoFilter(
        jPh     = jPh,
        jGt     = jGt,
        filtObj = traitsToKeep
    )
    tmpDF <- jTasFilt$phenoDf # check for missing values

    # Logic - Check for out of range p-values
    if (maxP > 1 || maxP < 0) {
        stop("p-value is out of range (0 - 1)")
    }

    # Logic - Convert p-values to Java data types
    if (!is.numeric(maxP)) {
        stop("p-value must be numeric")
    } else {
        maxP <- as.double(maxP)
    }

    # Logic - Convert threads to Java data types
    if (!is.null(maxThreads)) {
        maxThreads <- rJava::.jnew("java/lang/Integer", toString(maxThreads))
    } else {
        maxThreads <- rJava::.jnull()
    }

    # Logic - Check output data
    if (!is.null(outputFile)) {
        saveToFile <- TRUE
    } else {
        saveToFile <- FALSE
        outputFile <- rJava::.jnull()
        # outputFile <- "void"
    }

    # Logic - Handle association types and output
    if (!fitMarkers & is.null(kinship)) {
        if (!fastAssociation) {
            message("Association Analysis : BLUEs")
            assocOut <- jRC$association(
                rJava::.jnull(),
                rJava::.jnull(),
                jTasFilt$phenotype,
                rJava::.jnull(),
                as.integer(minClassSize),
                biallelicOnly,
                appendAddDom,
                saveToFile,
                outputFile,
                maxP
            )
            assocType <- "BLUE"
        } else {
            stop("Don't know how to analyze with given parameter inputs.")
        }
    } else if (fitMarkers & is.null(kinship)) {
        if (!fastAssociation) {
            if (!rJava::is.jnull(jGt)) {
                if (any(jTasFilt$attTypes == "factor")) {
                    message("Association Analysis : GLM")
                    message("(NOTE) Factors detected - running initial BLUE calculation...")
                    blueOut <- jRC$association(
                        rJava::.jnull(),
                        rJava::.jnull(),
                        jTasFilt$phenotype,
                        rJava::.jnull(),
                        as.integer(minClassSize),
                        biallelicOnly,
                        appendAddDom,
                        saveToFile,
                        outputFile,
                        maxP
                    )
                    blueOut <- blueOut$get("BLUE")
                    message("(NOTE) BLUEs calculated - using output to test markers...")
                    blueOut <- combineTasselGenotypePhenotype(
                        genotypeTable = jGt,
                        phenotype = blueOut
                    )
                    assocOut <- jRC$association(
                        rJava::.jnull(),
                        blueOut$genotypeTable(),
                        blueOut$phenotype(),
                        blueOut,
                        as.integer(minClassSize),
                        biallelicOnly,
                        appendAddDom,
                        saveToFile,
                        outputFile,
                        maxP
                    )
                    assocType <- "GLM"
                } else {
                    message("Association Analysis : GLM")
                    assocOut <- jRC$association(
                        rJava::.jnull(),
                        jTasFilt$genotypeTable,
                        jTasFilt$phenotype,
                        jTasFilt$genotypePhenotype,
                        as.integer(minClassSize),
                        biallelicOnly,
                        appendAddDom,
                        saveToFile,
                        outputFile,
                        maxP
                    )
                    assocType <- "GLM"
                }
            } else {
                stop("`assocModelFitter()` needs genotype data to fit markers")
            }
        } else {
            if (!rJava::is.jnull(jGt)) {
                if (any(apply(tmpDF, 2, function(x) any(is.na(x))))) {
                    stop("Missing phenotype data entries detected!")
                } else if (any(jTasFilt$attTypes == "factor")) {
                    message("Association Analysis : Fast Association")
                    message("(NOTE) Factors detected - running initial BLUE calculation...")
                    blueOut <- jRC$association(
                        rJava::.jnull(),
                        rJava::.jnull(),
                        jTasFilt$phenotype,
                        rJava::.jnull(),
                        as.integer(minClassSize),
                        biallelicOnly,
                        appendAddDom,
                        saveToFile,
                        outputFile,
                        maxP
                    )
                    blueOut <- blueOut$get("BLUE")
                    message("(NOTE) BLUEs calculated - using output to test markers...")
                    blueOut <- combineTasselGenotypePhenotype(
                        genotypeTable = jGt,
                        phenotype = blueOut
                    )

                    if (rJava::is.jnull(outputFile)) {
                        outputFile <- "void"
                        # saveToFile <- FALSE
                    }
                    assocOut <- rJava::.jcall(
                        "net.maizegenetics.plugindef.GenerateRCode", # fully‑qualified class
                        "Ljava/util/Map;",                           # JNI return type
                        "fastAssociation",                           # static method name
                        blueOut,                                     # GenotypePhenotype Java object
                        as.double(maxP),                             # primitive double (maxp)
                        maxThreads,                                  # max threads
                        saveToFile,                                  # primitive boolean (writeToFile)
                        outputFile                                   # java.lang.String (outputFile)
                    )
                    assocType <- "FastAssoc"
                } else {
                    if (rJava::is.jnull(outputFile)) {
                        outputFile <- "void"
                        # saveToFile <- FALSE
                    }
                    message("Association Analysis : Fast Association")
                    assocOut <- rJava::.jcall(
                        "net.maizegenetics.plugindef.GenerateRCode", # fully‑qualified class
                        "Ljava/util/Map;",                           # JNI return type
                        "fastAssociation",                           # static method name
                        jTasFilt$genotypePhenotype,                  # GenotypePhenotype Java object
                        as.double(maxP),                             # primitive double (maxp)
                        maxThreads,                                  # java.lang.Integer (maxThreads)
                        saveToFile,                                  # primitive boolean (writeToFile)
                        outputFile                                   # java.lang.String (outputFile)
                    )
                    assocType <- "FastAssoc"
                }
            } else {
                stop("`assocModelFitter()` needs genotype data to fit markers")
            }
        }
    } else if (fitMarkers & !is.null(kinship) & !fastAssociation) {
        if (!rJava::is.jnull(jGt)) {
            if (any(jTasFilt$attTypes == "factor")) {
                message("Association Analysis : MLM")
                message("(NOTE) Factors detected - running initial BLUE calculation...")
                blueOut <- jRC$association(
                    rJava::.jnull(),
                    rJava::.jnull(),
                    jTasFilt$phenotype,
                    rJava::.jnull(),
                    as.integer(minClassSize),
                    biallelicOnly,
                    appendAddDom,
                    saveToFile,
                    outputFile,
                    maxP
                )
                blueOut <- blueOut$get("BLUE")
                message("(NOTE) BLUEs calculated - using output to test markers...")
                blueOut <- combineTasselGenotypePhenotype(
                    genotypeTable = jGt,
                    phenotype = blueOut
                )
                assocOut <- jRC$association(
                    kinship,
                    blueOut$genotypeTable(),
                    blueOut$phenotype(),
                    blueOut,
                    as.integer(minClassSize),
                    biallelicOnly,
                    appendAddDom,
                    saveToFile,
                    outputFile,
                    maxP
                )
                assocType <- "MLM"
            } else {
                message("Association Analysis : MLM")
                assocOut <- jRC$association(
                    kinship,
                    jTasFilt$genotypeTable,
                    jTasFilt$phenotype,
                    jTasFilt$genotypePhenotype,
                    as.integer(minClassSize),
                    biallelicOnly,
                    appendAddDom,
                    saveToFile,
                    outputFile,
                    maxP
                )
                assocType <- "MLM"
            }
        } else {
            stop("`assocModelFitter()` needs genotype data to fit markers")
        }
    } else {
        stop("Don't know how to analyze with given parameter inputs.")
    }

    # DEBUG
    # return(tableReportList(assocOut))

    if (!saveToFile) {
        return(
            tableReportListToAssociationResults(
                trl   = tableReportList(assocOut),
                aType = assocType
            )
        )
    } else {
        message("Saved output to disk")
        return(NULL)
    }

}



## Phenotype filter - return modified TASSEL object - not exported (house keeping)
##
## 'jPh' is a Java Phenotype and 'jGt' a Java GenotypeTable (possibly a Java
## null), both as resolved by '.resolveTasselInput()'.
#' @importFrom rlang .data
tasPhenoFilter <- function(jPh, jGt, filtObj) {

    # Get all TASSEL object trait metadata
    phenoAttDf <- extractPhenotypeAttDf(jPh)

    # Get phenotype data frame
    phenoDF <- tableReportToDF(jPh)

    # Convert <data> and <covariates> to doubles (correct pass to TASSEL)
    doubCols <- as.character(
        phenoAttDf$traitName[which(phenoAttDf$traitType == "data" | phenoAttDf$traitType == "covariate")]
    )
    phenoDF[doubCols] <- sapply(phenoDF[doubCols], as.double)

    # Get taxa column
    taxaCol <- as.character(phenoAttDf$traitName[which(phenoAttDf$traitType == "taxa")])
    taxaNames <- as.vector(phenoDF[[taxaCol]])

    # Get non-taxa columns and reorder filtered columns (correct pass to TASSEL)
    origColNames <- colnames(phenoDF)
    filtObjRight <- c(taxaCol, filtObj)
    filtObjRight <- filtObjRight[match(origColNames, filtObjRight)]
    filtObjRight <- filtObjRight[!is.na(filtObjRight)]

    # Filter data frame columns based on association formula
    phenoDF <- phenoDF[, filtObjRight]
    phenoAttDf <- phenoAttDf[phenoAttDf$traitName %in% filtObjRight, , drop = FALSE]

    # Get vector of non-taxa column names
    phenoColNames <- colnames(phenoDF)
    notTaxaCols <- phenoColNames[!(phenoColNames %in% taxaCol)]

    # Get attribute types
    attTypes <- as.vector(phenoAttDf$traitType[which(phenoAttDf$traitType != "taxa")])

    # Send filtered data frame to TASSEL methods
    jList <- rJava::new(rJava::J("java/util/ArrayList"))
    for (col_i in notTaxaCols) {
        jList$add(rJava::.jarray(phenoDF[[col_i]]))
    }

    # Create filtered TASSEL phenotype
    jc <- rJava::J("net/maizegenetics/plugindef/GenerateRCode")
    jc <- jc$createPhenotypeFromRDataFrameElements(
        taxaNames,
        rJava::.jarray(notTaxaCols),
        rJava::.jarray(attTypes),
        jList
    )

    # Return modified TASSEL objects
    if (rJava::is.jnull(jGt)) {
        return(
            list(
                attTypes = attTypes,
                phenoDf = phenoDF,
                finalVars = filtObj,
                notTaxaCols = notTaxaCols,
                genotypeTable = rJava::.jnull(),
                phenotype = jc,
                genotypePhenotype = rJava::.jnull()
            )
        )
    } else {
        jcComb <- combineTasselGenotypePhenotype(
            genotypeTable = jGt,
            phenotype = jc
        )
        return(
            list(
                attTypes = attTypes,
                phenoDf = phenoDF,
                finalVars = filtObj,
                notTaxaCols = notTaxaCols,
                genotypeTable = jcComb$genotypeTable(),
                phenotype = jcComb$phenotype(),
                genotypePhenotype = jcComb
            )
        )
    }
}


