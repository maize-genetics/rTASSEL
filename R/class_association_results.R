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
#'
#' @name AssociationResults-class
#' @rdname AssociationResults-class
#' @exportClass AssociationResults
setClass(
    Class = "AssociationResults",
    representation = representation(
        results = "list",
        traits = "character"
    )
)

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


