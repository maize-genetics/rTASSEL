## ----
#' @title LDRegion Class
#'
#' @description An S4 class to define a genomic region for highlighting
#'   LD blocks on a \code{\link{plotLD}} plot. Each region is specified
#'   by a start and end position (in base pairs) with an optional text
#'   label, outline color, line thickness, and span display control.
#'
#' @slot start A single numeric value for the region start position (bp).
#' @slot end A single numeric value for the region end position (bp).
#'   Must be greater than or equal to \code{start}.
#' @slot label A single character string used as the block annotation.
#'   Defaults to \code{NA_character_} (no label).
#' @slot color A single character string specifying the outline color
#'   for the block highlight. Must be a valid R color. Defaults to
#'   \code{"black"}.
#' @slot linewidth A single numeric value controlling the outline
#'   thickness (in mm) for the block highlight. Defaults to
#'   \code{NA_real_}, which uses the plot-level auto-calculated width.
#' @slot showSpan Logical. If \code{TRUE}, the genomic span (e.g.,
#'   \dQuote{342.9 kbp}) is shown in the block annotation regardless
#'   of whether a \code{label} is provided. If \code{FALSE}, the span
#'   is never shown. Defaults to \code{TRUE}.
#'
#' @seealso \code{\link{plotLD}}
#'
#' @name LDRegion-class
#' @rdname LDRegion-class
#' @exportClass LDRegion
setClass(
    Class = "LDRegion",
    slots = c(
        start     = "numeric",
        end       = "numeric",
        label     = "character",
        color     = "character",
        linewidth = "numeric",
        showSpan  = "logical"
    ),
    prototype = list(
        start     = NA_real_,
        end       = NA_real_,
        label     = NA_character_,
        color     = "black",
        linewidth = NA_real_,
        showSpan  = TRUE
    )
)


## ----
#' @title LDRegion validation
#'
#' @name LDRegion-validity
#'
#' @description
#' Validates that an \code{LDRegion} object has properly formed slots:
#' \code{start} and \code{end} must be length-1 finite numerics with
#' \code{end >= start}, \code{color} must be a recognized R color,
#' \code{linewidth} must be length-1 (positive or NA), and
#' \code{showSpan} must be a length-1 non-NA logical.
#'
#' @param object An \code{LDRegion} object.
setValidity("LDRegion", function(object) {
    errors <- character()

    if (length(object@start) != 1 || !is.finite(object@start)) {
        errors <- c(errors, "'start' must be a single finite numeric value")
    }
    if (length(object@end) != 1 || !is.finite(object@end)) {
        errors <- c(errors, "'end' must be a single finite numeric value")
    }
    if (length(errors) == 0 && object@end < object@start) {
        errors <- c(errors, "'end' must be greater than or equal to 'start'")
    }

    if (length(object@label) != 1) {
        errors <- c(errors, "'label' must be a single character value or NA")
    }

    if (length(object@color) != 1 || is.na(object@color)) {
        errors <- c(errors, "'color' must be a single non-NA character value")
    } else {
        validColor <- tryCatch(
            { grDevices::col2rgb(object@color); TRUE },
            error = function(e) FALSE
        )
        if (!validColor) {
            errors <- c(errors, paste0(
                "'color' (\"", object@color, "\") is not a valid R color"
            ))
        }
    }

    if (length(object@linewidth) != 1) {
        errors <- c(errors, "'linewidth' must be a single numeric value or NA")
    } else if (!is.na(object@linewidth) && object@linewidth <= 0) {
        errors <- c(errors, "'linewidth' must be positive")
    }

    if (length(object@showSpan) != 1 || is.na(object@showSpan)) {
        errors <- c(errors, "'showSpan' must be TRUE or FALSE")
    }

    if (length(errors) == 0) TRUE else errors
})


## ----
#' @title Create an LDRegion object
#'
#' @description Constructor for \code{\linkS4class{LDRegion}} objects.
#'   Defines a single genomic region for LD block highlighting.
#'
#' @param start Numeric. Start position of the region in base pairs.
#' @param end Numeric. End position of the region in base pairs.
#'   Must be \code{>= start}.
#' @param label Character. Optional text label for the block.
#'   Defaults to \code{NA_character_} (no label).
#' @param color Character. Outline color for the block highlight.
#'   Must be a valid R color. Defaults to \code{"black"}.
#' @param linewidth Numeric. Outline thickness (mm) for the block
#'   highlight. \code{NA} (the default) uses the plot-level
#'   auto-calculated width.
#' @param showSpan Logical. If \code{TRUE} (default), the genomic span
#'   (e.g., \dQuote{342.9 kbp}) is appended to the block annotation.
#'   When no \code{label} is provided, the span is shown on its own.
#'   Set to \code{FALSE} to suppress the span entirely.
#'
#' @return An object of class \code{\linkS4class{LDRegion}}.
#'
#' @examples
#' \dontrun{
#' # Region with a label (span shown by default)
#' LDRegion(start = 157104, end = 500000, label = "Block A")
#'
#' # Region with custom color, no label but span still shown
#' LDRegion(start = 800000, end = 1200000, color = "blue")
#'
#' # Suppress the span text
#' LDRegion(start = 100, end = 500, label = "QTL", showSpan = FALSE)
#'
#' # Custom line thickness
#' LDRegion(start = 100, end = 500, linewidth = 1.5)
#'
#' # Pass a list of regions to plotLD
#' plotLD(
#'   tasObj,
#'   ldBlocks = list(
#'     LDRegion(start = 157104, end = 500000, label = "A"),
#'     LDRegion(start = 800000, end = 1200000, label = "B", color = "blue")
#'   )
#' )
#' }
#'
#' @importFrom methods new validObject
#' @importFrom grDevices col2rgb
#' @export
LDRegion <- function(
    start,
    end,
    label     = NA_character_,
    color     = "black",
    linewidth = NA_real_,
    showSpan  = TRUE
) {
    methods::new(
        "LDRegion",
        start     = as.numeric(start),
        end       = as.numeric(end),
        label     = as.character(label),
        color     = as.character(color),
        linewidth = as.numeric(linewidth),
        showSpan  = as.logical(showSpan)
    )
}


# /// Methods (show) ////////////////////////////////////////////////

## ----
#' @title Show method for LDRegion objects
#'
#' @description Prints a compact summary of an \code{LDRegion} object.
#'
#' @param object An \code{\linkS4class{LDRegion}} object.
#'
#' @rdname LDRegion-class
#' @aliases show,LDRegion-method
setMethod(
    f = "show",
    signature = "LDRegion",
    definition = function(object) {
        pointerSymbol <- cli::col_green(cli::symbol$pointer)

        spanBp <- object@end - object@start

        spanStr <- if (spanBp >= 1e6) {
            paste0(round(spanBp / 1e6, 2), " Mbp")
        } else if (spanBp >= 1e3) {
            paste0(round(spanBp / 1e3, 2), " kbp")
        } else {
            paste0(spanBp, " bp")
        }

        labelStr <- if (!is.na(object@label)) {
            cli::style_bold(object@label)
        } else {
            cli::col_grey("NA")
        }

        lwStr <- if (is.na(object@linewidth)) {
            cli::col_grey("auto")
        } else {
            paste0(object@linewidth, " mm")
        }

        fmtBp <- function(x) format(x, big.mark = ",", scientific = FALSE)
        rangeStr <- paste0(
            cli::style_bold(fmtBp(object@start)),
            " - ",
            cli::style_bold(fmtBp(object@end)),
            " bp (",
            cli::style_bold(spanStr),
            ")"
        )

        msg <- c(
            paste0("An ", cli::style_bold("LDRegion"), " object"),
            paste0(" ", pointerSymbol, " Range......: ", rangeStr),
            "---",
            paste0(" ", pointerSymbol, " Label......: ", labelStr),
            paste0(" ", pointerSymbol, " Color......: ", cli::style_bold(object@color)),
            paste0(" ", pointerSymbol, " Linewidth..: ", lwStr),
            paste0(" ", pointerSymbol, " Show span..: ", cli::style_bold(object@showSpan))
        )

        cat(msg, sep = "\n")
    }
)
