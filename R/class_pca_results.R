## ----
#' @title PCAResults Class
#'
#' @description
#' Class \code{PCAResults} defines a \code{rTASSEL}
#' Class for storing TASSEL 5 PCA results
#'
#' @slot results A list of \code{data.frame} objects containing summary results
#' @slot jObj An rJava reference object pointing to PCA results in Java memory
#'
#' @name PCAResults-class
#' @rdname PCAResults-class
#' @exportClass PCAResults
setClass(
    Class = "PCAResults",
    representation = representation(
        results = "list",
        jObj = "jobjRef"
    )
)


## ----
#' @title PCAResults validation
#'
#' @name PCAResults-validity
#'
#' @description
#' Checks if all elements of list are \code{data.frame} objects
#'
#' @param object A \code{PCAResults} object
setValidity("PCAResults", function(object) {
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
#' @title Show methods for PCAResults objects
#'
#' @description
#' Prints outs summary information from \code{PCAResults} objects
#'
#' @param object a \code{\linkS4class{PCAResults}} object
#'
#' @docType methods
#' @rdname PCAResults-class
#' @aliases show,PCAResults-method
setMethod(
    f = "show",
    signature = "PCAResults",
    definition = function(object) {
        # local parameters
        maxTraitsToPrint <- 5
        indentStyle <- "  *"

        # header text
        titleMsg <- paste(
            "PCAResults object with", length(object@results),
            "reports and", ncol(object@results$PC_Datum) - 1, "reported PCs \n"
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
    signature = "PCAResults",
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
        assocRes   = "PCAResults",
        reportName = "ANY"
    ),
    definition = function(assocRes, reportName) {
        if (missing(reportName)) reportName <- NULL
        returnReportElements(
            assocRes             = assocRes,
            reportName           = reportName,
            defaultReportElement = "PC_Datum"
        )
    }
)


