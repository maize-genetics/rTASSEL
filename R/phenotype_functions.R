## ----
#' @title Wrapper function of TasselGenotypePhenotype class for phenotype
#'    data from a path.
#'
#' @description
#' \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
#'
#' This function is a wrapper for the deprecated
#' \code{TasselGenotypePhenotype} class. It is used for storing phenotype
#' information into a class object. This will read in phenotype data from
#' a path. Use \code{\link{readPhenotype}()} instead, which returns a
#' \code{\linkS4class{TasselPhenotype}} object.
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @name readPhenotypeFromPath
#' @rdname readPhenotypeFromPath
#'
#' @param path A phenotype data path.
#'
#' @seealso \code{\link{readPhenotype}}
#'
#' @importFrom rJava J
#' @importFrom rJava %instanceof%
#' @importFrom rJava new
#' @export
readPhenotypeFromPath <- function(path) {
    lifecycle::deprecate_warn(
        "0.14.0", "readPhenotypeFromPath()", "readPhenotype()"
    )

    if (!file.exists(path)) {
        stop("Cannot open file ", path, ": No such file or directory")
    }

    jObj <- rJava::new(
        rJava::J("net.maizegenetics.phenotype.PhenotypeBuilder")
    )$fromFile(path)

    return(.tasselObjectConstructor(jObj$build()$get(0L)))
}


#' @title Wrapper function of TasselGenotypePhenotype class for phenotype
#'    data from an R data frame
#'
#' @description
#' \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
#'
#' This function is a wrapper for the deprecated
#' \code{TasselGenotypePhenotype} class. It is used for storing phenotype
#' information into a class object. Use \code{\link{readPhenotype}()}
#' instead, which returns a \code{\linkS4class{TasselPhenotype}} object.
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @name readPhenotypeFromDataFrame
#' @rdname readPhenotypeFromDataFrame
#'
#' @param phenotypeDF A \code{R} object of class \code{data.frame}.
#' @param taxaID The column name that represents your taxa data as a string.
#' @param attributeTypes A vector of non-taxa attributes. If \code{NULL}, all
#'    attributes will be TASSEL \code{<data>} types.
#'
#' @seealso \code{\link{readPhenotype}}
#'
#' @importFrom rJava .jarray
#' @importFrom rJava J
#' @importFrom rJava new
#' @export
readPhenotypeFromDataFrame <- function(phenotypeDF,
                                       taxaID,
                                       attributeTypes = NULL) {
    lifecycle::deprecate_warn(
        "0.14.0", "readPhenotypeFromDataFrame()", "readPhenotype()"
    )

    safeAtt <- c("covariate", "data", "factor", "taxa")
    if (!is.null(attributeTypes) & !all(attributeTypes %in% safeAtt)) {
        stop(
            paste0(
                "Parameter `attributeTypes` contains incorrect attributes.\n",
                "Please select from the following:\n",
                "  taxa\n",
                "  factor\n",
                "  data\n",
                "  covariate\n"
            )
        )
    }

    # TODO Remove tibble check
    if (inherits(phenotypeDF, "tbl_df")) {
        phenotypeDF <- as.data.frame(phenotypeDF)
    }

    taxaNames <- as.vector(phenotypeDF[, taxaID])
    colnames <- colnames(phenotypeDF)
    notTaxaCols <- colnames[!colnames %in% taxaID]
    if(is.null(attributeTypes)) {
        atttype <- c(rep("data", length(notTaxaCols)))
    } else {
        atttype <- attributeTypes
    }
    jList <- rJava::new(rJava::J("java/util/ArrayList"))
    for (col_i in notTaxaCols) {
        jList$add(.jarray(phenotypeDF[[col_i]]))
    }
    jc <- J("net/maizegenetics/plugindef/GenerateRCode")
    jc <- jc$createPhenotypeFromRDataFrameElements(
        taxaNames,
        rJava::.jarray(notTaxaCols),
        rJava::.jarray(atttype),
        jList
    )
    return(.tasselObjectConstructor(jc))
}


#' @title Get an R/\code{DataFrame} phenotype data frame from TASSEL object
#'
#' @description
#' \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
#'
#' This function will extract a \code{DataFrame}-based R data
#'    frame from an object that contains phenotypic data. Use
#'    \code{\link[base]{as.data.frame}()} on a
#'    \code{\linkS4class{TasselPhenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}} instead. Column data will be
#'    converted to the following
#'    types data depending on TASSEL data type:
#'    \itemize{
#'      \item{\code{<taxa>}: \code{character}}
#'      \item{\code{<data>}: \code{numeric}}
#'      \item{\code{<covariate>}: \code{numeric}}
#'      \item{\code{<factor>}: \code{factor}}
#'    }
#'
#' @return Returns an \code{DataFrame} based data frame
#'
#' @name getPhenotypeDF
#' @rdname getPhenotypeDF
#'
#' @param tasObj An object of class \code{\linkS4class{TasselPhenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the deprecated
#'    \code{TasselGenotypePhenotype} class are still accepted.
#'
#' @seealso \code{\link[base]{as.data.frame}}
#'
#' @importFrom rJava is.jnull
#' @export
getPhenotypeDF <- function(tasObj) {
    lifecycle::deprecate_warn(
        "0.14.0", "getPhenotypeDF()", I("`as.data.frame()`")
    )

    jPhenoTable <- .resolveTasselInput(
        tasObj, "phenotype", "getPhenotypeDF"
    )$jPh

    jPhenoAttri <- extractPhenotypeAttDf(jPhenoTable)
    jPhenoTable <- tableReportToDF(jPhenoTable)

    # # Get list of TASSEL data types
    # attributes <- c("taxa", "factor", "covariate", "data")
    # att <- lapply(seq_along(attributes), function(i) {
    #     as.vector(jPhenoAttri$traitName[which(jPhenoAttri$traitType == attributes[i])])
    # })
    # names(att) <- attributes
    #
    # # Convert column data based TASSEL data types
    # jPhenoTable[c(att$covariate, att$data)] <- sapply(
    #     jPhenoTable[c(att$covariate, att$data)], as.numeric
    # )
    # jPhenoTable[c(att$factor)] <- lapply(
    #     jPhenoTable[c(att$factor)], factor
    # )
    # names(jPhenoTable) <- jPhenoAttri$traitName
    return(jPhenoTable)
}


## Get a Phenotype object - not exported (house keeping)
getPhenotypeTable <- function(jtsObject) {
    if(is(jtsObject, "TasselGenotypePhenotype")) {
        return(jtsObject@jPhenotypeTable)
    }
    jtsObject <- .unwrapTasselObject(jtsObject)
    if(!is(jtsObject,"jobjRef")) return(rJava::.jnull())
    if(jtsObject %instanceof% "net.maizegenetics.phenotype.Phenotype") {
        return(jtsObject)
    } else if(jtsObject %instanceof% "net.maizegenetics.phenotype.GenotypePhenotype") {
        return(jtsObject$phenotype())
    } else {
        return(rJava::.jnull())
    }
}


## Get Phenotype attributes as data frame - not exported (house keeping)
extractPhenotypeAttDf <- function(phenotype) {
    traitName = phenotype$getTableColumnNames()
    traitType = unlist(
        lapply(as.list(phenotype$typeListCopy()), function(tc) {
            tc$toString()
        })
    )

    # Pull the java class and return the class without the whole path
    traitAttribute = unlist(
        lapply(as.list(phenotype$attributeListCopy()), function(tc) {
            strsplit(tc$getClass()$toString(),"\\.")[[1]][4]
        })
    )
    return(tibble::tibble(traitName, traitType, traitAttribute))
}
