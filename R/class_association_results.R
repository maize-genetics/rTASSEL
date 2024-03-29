## ----
#' @title AssociationResults Class
#'
#' @description
#' Class \code{AssociationResults} defines a \code{rTASSEL}
#' Class for storing TASSEL 5 and general-purpose association results
#' from *WAS studies.
#'
#' @slot results A list of \code{data.frame} objects containing summary results
#' @slot traits A vector of type \code{character} containing the trait IDs
#'    modeled.
#' @slot assocType \code{character} object describing association type
#'
#' @name AssociationResults-class
#' @rdname AssociationResults-class
#' @exportClass AssociationResults
setClass(
    Class = "AssociationResults",
    representation = representation(
        results = "list",
        traits = "character",
        assocType = "character"
    )
)


## ----
#' @title AssociationResults validation
#'
#' @name AssociationResults-validity
#'
#' @description
#' Checks if all elements of list are \code{data.frame} objects
#'
#' @param object A \code{AssociationResults} object
setValidity("AssociationResults", function(object) {
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
#' @title Show methods for AssociationResults objects
#'
#' @description
#' Prints outs summary information from \code{AssociationResults} objects
#'
#' @param object a \code{\linkS4class{AssociationResults}} object
#'
#' @docType methods
#' @rdname AssociationResults-class
#' @aliases show,AssociationResults-method
setMethod(
    f = "show",
    signature = "AssociationResults",
    definition = function(object) {
        # local parameters
        maxTraitsToPrint <- 5
        indentStyle <- "  *"

        # header text
        titleMsg <- paste(
            "AssociationResults object with", length(object@results),
            "reports and", length(object@traits), "mapped traits \n"
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

        # trait text
        cat("Traits:\n")
        if (length(object@traits) <= maxTraitsToPrint) {
            for (i in seq_len(length(object@traits))) {
                cat(indentStyle, object@traits[i], "\n")
            }
        } else {
            remTraits <- length(object@traits) - maxTraitsToPrint

            for (i in seq_len(maxTraitsToPrint)) {
                cat(indentStyle, object@traits[i], "\n")
            }
            cat("    ...\n")
            cat(paste("  (with", remTraits, "more values)\n"))
        }
    }
)


## ----
#' @rdname reportNames
#' @export
setMethod(
    f = "reportNames",
    signature = "AssociationResults",
    definition = function(object) {
        return(names(object@results))
    }
)


## ----
#' @rdname reportNames
#' @export
setMethod(
    f = "traitNames",
    signature = "AssociationResults",
    definition = function(object) {
        return(object@traits)
    }
)


## ----
#' @rdname tableReport
#' @export
setMethod(
    f = "tableReport",
    signature = signature(
        assocRes   = "AssociationResults",
        reportName = "ANY"
    ),
    definition = function(assocRes, reportName) {
        if (missing(reportName)) {
            reportName <- NULL
        }

        if (!is.character(reportName) && !is.null(reportName)) {
            stop("'reportName' must be of type 'character'")
        }

        if (is.null(reportName)) {
            return(assocRes@results)
        } else {
            if (!reportName %in% reportNames(assocRes)) {
                stop("Report ID not found in object")
            }
            return(assocRes@results[[reportName]])
        }
    }
)


## ----
#' @rdname associationType
#' @export
setMethod(
    f = "associationType",
    signature = signature("AssociationResults"),
    definition = function(assocRes) {
        return(assocRes@assocType)
    }
)


