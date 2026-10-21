## ----
#' @title MDSResults Class
#'
#' @description
#' Class \code{MDSResults} defines a \code{rTASSEL}
#' Class for storing TASSEL 5 MDS results
#'
#' @slot results A named list of \code{tibble} objects containing summary
#' results
#' @slot jObj An rJava reference object pointing to MDS results in Java memory
#'
#' @name MDSResults-class
#' @rdname MDSResults-class
#' @exportClass MDSResults
setClass(
    Class = "MDSResults",
    representation = representation(
        results = "list",
        jObj = "jobjRef"
    )
)


## ----
#' @title MDSResults validation
#'
#' @name MDSResults-validity
#'
#' @description
#' Checks if all elements of list are \code{data.frame} objects
#'
#' @param object A \code{MDSResults} object
setValidity("MDSResults", function(object) {
    errors <- character()

    dfCheck <- vapply(object@results, is.data.frame, logical(1))

    if (!all(dfCheck)) {
        invalidIndices <- which(!dfCheck)
        if (length(invalidIndices) == 1) {
            grammar <- c("Element", "is")
        } else {
            invalidIndices <- paste(invalidIndices, collapse = ", ")
            grammar <- c("Elements", "are")
        }
        invalidElements <- paste(
            grammar[1], invalidIndices, grammar[2], "not a data frame."
        )
        msg <- paste0("Invalid list: ", invalidElements)
        errors <- c(errors, msg)
    }

    if (length(errors) == 0) {
        return(TRUE)
    } else {
        return(errors)
    }
})


## ----
#' @title Show methods for MDSResults objects
#'
#' @description
#' Prints outs summary information from \code{MDSResults} objects
#'
#' @param object a \code{\linkS4class{MDSResults}} object
#'
#' @docType methods
#' @rdname MDSResults-class
#' @aliases show,MDSResults-method
setMethod(
    f = "show",
    signature = "MDSResults",
    definition = function(object) {
        # local parameters
        indentStyle <- "  *"

        # header text
        titleMsg <- paste(
            "MDSResults object with", length(object@results),
            "reports and", ncol(object@results$MDS_PCs_Datum) - 1,
            "reported axes \n"
        )
        cat(titleMsg)

        # results text
        cat("Results:\n")
        dims <- vapply(object@results, dim, numeric(2))
        for (i in seq_along(object@results)) {
            dfName  <- names(object@results[i])
            currDim <- paste0("(", paste(dims[, dfName], collapse = ", "), ")")
            cat(indentStyle, dfName, currDim, "\n")
        }
    }
)


## ----
#' @rdname reportNames
#' @export
setMethod(
    f = "reportNames",
    signature = "MDSResults",
    definition = function(object) {
        return(names(object@results))
    }
)


## ----
#' @rdname tableReport
#' @export
setMethod(
    f = "tableReport",
    signature = signature(
        assocRes   = "MDSResults",
        reportName = "ANY"
    ),
    definition = function(assocRes, reportName) {
        if (missing(reportName)) reportName <- NULL
        returnReportElements(
            results              = assocRes@results,
            reportName           = reportName,
            defaultReportElement = "MDS_PCs_Datum"
        )
    }
)


