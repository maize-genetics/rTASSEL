## ----
#' @title LDResults Class
#'
#' @description An S4 class that stores the results of a linkage
#'   disequilibrium (LD) analysis produced by \code{\link{linkageDiseq}}.
#'   The object contains the full pairwise LD table together with the
#'   parameters used to generate it.
#'
#' @slot results A \code{data.frame} containing the pairwise LD
#'   statistics. Expected columns include \code{Locus1}, \code{Position1},
#'   \code{Locus2}, \code{Position2}, \code{R^2}, \code{DPrime},
#'   \code{pDiseq}, and \code{N}, among others.
#' @slot ldType A single character string indicating the LD calculation
#'   type: \code{"All"} or \code{"SlidingWindow"}.
#' @slot windowSize A single numeric value for the sliding-window size.
#'   \code{NA_real_} when \code{ldType} is \code{"All"}.
#' @slot hetCalls A single character string indicating how heterozygous
#'   calls were handled: \code{"missing"}, \code{"ignore"}, or
#'   \code{"third"}.
#'
#' @seealso \code{\link{linkageDiseq}}, \code{\link{plotLD}}
#'
#' @name LDResults-class
#' @rdname LDResults-class
#' @exportClass LDResults
setClass(
    Class = "LDResults",
    slots = c(
        results    = "data.frame",
        ldType     = "character",
        windowSize = "numeric",
        hetCalls   = "character"
    )
)


## ----
#' @title LDResults validation
#'
#' @name LDResults-validity
#'
#' @description
#' Validates that an \code{LDResults} object has a properly formed
#' \code{results} data frame with the required LD columns and that
#' \code{ldType} and \code{hetCalls} contain allowed values.
#'
#' @param object An \code{LDResults} object.
setValidity("LDResults", function(object) {
    errors <- character()

    requiredCols <- c(
        "Locus1", "Position1", "Locus2", "Position2",
        "R^2", "DPrime", "pDiseq", "N"
    )
    missing <- setdiff(requiredCols, colnames(object@results))
    if (length(missing) > 0) {
        errors <- c(errors, paste0(
            "'results' is missing required columns: ",
            paste(missing, collapse = ", ")
        ))
    }

    if (length(object@ldType) != 1 || !object@ldType %in% c("All", "SlidingWindow")) {
        errors <- c(errors, "'ldType' must be \"All\" or \"SlidingWindow\"")
    }

    if (length(object@windowSize) != 1) {
        errors <- c(errors, "'windowSize' must be a single numeric value or NA")
    }

    if (length(object@hetCalls) != 1 || !object@hetCalls %in% c("missing", "ignore", "third")) {
        errors <- c(errors, "'hetCalls' must be \"missing\", \"ignore\", or \"third\"")
    }

    if (length(errors) == 0) TRUE else errors
})


# /// Methods (show) ////////////////////////////////////////////////

## ----
#' @title Show method for LDResults objects
#'
#' @description Prints a compact summary of an \code{LDResults} object.
#'
#' @param object An \code{\linkS4class{LDResults}} object.
#'
#' @rdname LDResults-class
#' @aliases show,LDResults-method
setMethod(
    f = "show",
    signature = "LDResults",
    definition = function(object) {
        pointerSymbol <- cli::col_green(cli::symbol$pointer)

        df   <- object@results
        dims <- paste0(
            format(nrow(df), big.mark = ","),
            " pairs x ",
            ncol(df), " columns"
        )

        loci <- unique(c(df$Locus1, df$Locus2))
        nLoci <- length(loci)

        wsStr <- if (is.na(object@windowSize)) {
            cli::col_grey("NA")
        } else if (object@windowSize == -1) {
            cli::col_grey("NA (all pairs)")
        } else {
            as.character(object@windowSize)
        }

        msg <- c(
            paste0("An ", cli::style_bold("LDResults"), " object"),
            paste0(" ", pointerSymbol, " Dimensions.: ", cli::style_bold(dims)),
            paste0(" ", pointerSymbol, " Chromosomes: ", cli::style_bold(nLoci)),
            "---",
            paste0(" ", pointerSymbol, " LD type....: ", cli::style_bold(object@ldType)),
            paste0(" ", pointerSymbol, " Window size: ", wsStr),
            paste0(" ", pointerSymbol, " Het. calls.: ", cli::style_bold(object@hetCalls))
        )

        cat(msg, sep = "\n")
    }
)


# /// Methods (tableReport) /////////////////////////////////////////

## ----
#' @rdname tableReport
#' @export
setMethod(
    f = "tableReport",
    signature = signature(
        assocRes   = "LDResults",
        reportName = "ANY"
    ),
    definition = function(assocRes, reportName) {
        tibble::as_tibble(assocRes@results)
    }
)
