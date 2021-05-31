#---------------------------------------------------------------------
# Script Name:   PhenotypeFunctions.R
# Description:   Functions to create TASSEL Phenotype
# Author:        Brandon Monier & Ed Buckler
# Created:       2018-11-26 at 11:14:36
# Last Modified: 2021-05-31 at 11:51:11
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions
#    necessary for reading in phenotype datasets into R and
#    extracting data from TASSEL phenotype objects.
#--------------------------------------------------------------------

#' @title Wrapper function of TasselGenotypePhenotype class for phenotype
#'    data from a path.
#'
#' @description This function is a wrapper for the
#'    \code{TasselGenotypePhenotype} class. It is used for storing phenotype
#'    information into a class object. This will read in phenotype data from
#'    a path.
#'
#' @return Returns an object of \code{TasselGenotypePhenotype} class.
#'
#' @name readPhenotypeFromPath
#' @rdname readPhenotypeFromPath
#'
#' @param path A phenotype data path.
#'
#' @importFrom rJava J
#' @importFrom rJava %instanceof%
#' @importFrom rJava new
#'
#' @export
readPhenotypeFromPath <- function(path) {
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
#' @description This function is a wrapper for the
#'    \code{TasselGenotypePhenotype} class. It is used for storing phenotype
#'    information into a class object.
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
#' @importFrom rJava .jarray
#' @importFrom rJava J
#' @importFrom rJava new
#'
#' @export
readPhenotypeFromDataFrame <- function(phenotypeDF,
                                       taxaID,
                                       attributeTypes = NULL) {
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
    if (any(class(phenotypeDF) == "tbl_df")) {
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


#' @title Read phenotype information from a BrAPI web service
#'
#' @description Reads phenotype information from BreedBase using a BrAPI rest
#'    client.
#'
#' @return Returns an object of \code{data.frame} class.
#'
#' @name readPhenotypeFromBrapi
#' @rdname readPhenotypeFromBrapi
#'
#' @param brapiObj An BrAPI connection object of \code{BrapiCon} class.
#' @param studyName Return a specified study. If \code{NULL}, all data will be
#'    returned. Defaults to \code{NULL}.
#' @param verbose Show messages to console? Defaults to \code{TRUE}.
#'
#' @importFrom rPHG brapiURL
#' @importFrom httr GET
#' @importFrom httr content
#' @importFrom jsonlite fromJSON
#'
#' @export
readPhenotypeFromBrapi <- function(brapiObj, studyName = NULL, verbose = TRUE) {
    if (verbose) message("Downloading data...")
    url <- rPHG::brapiURL(brapiObj)

    obs <- paste0(url, "/observations?pageSize=10000")
    studies <- paste0(url, "/studies")

    ## Return observations
    xObs <- httr::GET(obs)
    xObs <- httr::content(xObs, as = "text", encoding = "ISO-8859-1")
    xObs <- jsonlite::fromJSON(xObs, flatten = TRUE)
    xObs <- xObs$result$data

    xObsSeason <- xObs$season
    xObsSeason <- do.call("rbind", xObsSeason)
    xObs <- xObs[, which(names(xObs) != "season")]
    xObs <- cbind(xObs, xObsSeason)


    ## Return studies
    xStu <- httr::GET(studies)
    xStu <- httr::content(xStu, as = "text", encoding = "ISO-8859-1")
    xStu <- jsonlite::fromJSON(xStu, flatten = TRUE)
    xStu <- xStu$result$data
    xStu <- xStu[!(colnames(xStu) %in% c("dataLinks", "seasons"))]

    ## JOIN study and observation info
    total <- merge(xObs, xStu, by = "studyDbId")

    ## Filter by parameters (WIP)
    if (!is.null(studyName)) total <- total[total$studyName %in% studyName, ]

    if (dim(total)[1] == 0) {
        warning("No data returned.", call. = FALSE)
        return(NULL)
    } else {
        return(total)
    }
}


#' @title Get an R/\code{DataFrame} phenotype data frame from TASSEL object
#'
#' @description This function will extract a \code{DataFrame}-based R data
#'    frame from an object of \code{TasselGenotypePhenotype} class that
#'    contains phenotypic data. Column data will be converted to the following
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
#' @param tasObj An object of class \code{TasselGenotypePenotype}.
#'
#' @importFrom rJava is.jnull
#'
#' @export
getPhenotypeDF <- function(tasObj) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jPhenoTable <- getPhenotypeTable(tasObj)
    if (rJava::is.jnull(jPhenoTable)) {
        stop("TASSEL phenotype object not found")
    }

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
    return(data.frame(traitName, traitType, traitAttribute))
}


