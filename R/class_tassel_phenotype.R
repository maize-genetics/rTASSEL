# /// S3 - method extensions ////////////////////////////////////////

methods::setOldClass("java_pheno_tbl")

## ----
#' @importFrom pillar pillar_shaft
#' @method pillar_shaft cov
#' @export
pillar_shaft.cov  <- function(x, ...) pillar::pillar_shaft(vctrs::vec_data(x), ...)

#' @importFrom pillar pillar_shaft
#' @method pillar_shaft data
#' @export
pillar_shaft.data <- function(x, ...) pillar::pillar_shaft(vctrs::vec_data(x), ...)

#' @importFrom pillar pillar_shaft
#' @method pillar_shaft fact
#' @export
pillar_shaft.fact <- function(x, ...) pillar::new_pillar_shaft_simple(paste0(cli::style_italic(x)), align = "left")

#' @importFrom pillar pillar_shaft
#' @method pillar_shaft taxa
#' @export
pillar_shaft.taxa <- function(x, ...) pillar::new_pillar_shaft_simple(paste0(cli::style_bold(x)), align = "left")


## ----
#' @importFrom pillar tbl_format_header
#' @method tbl_format_header java_pheno_tbl
#' @export
tbl_format_header.java_pheno_tbl <- function(x, ...) {
    nTaxa <- attr(x, "nTaxa")
    nTraits <- attr(x, "nTraits")
    header <- sprintf(
        "# A %s object: %s taxa %s %s traits",
        cli::style_bold("TasselPhenotype"),
        nTaxa,
        cli::symbol$times,
        nTraits
    )

    return(pillar::style_subtle(header))
}


## ----
#' @importFrom pillar tbl_format_footer
#' @method tbl_format_footer java_pheno_tbl
#' @export
tbl_format_footer.java_pheno_tbl <- function(x, ...) {
    defaultFooter <- NextMethod()
    nDfRow <- attr(x, "nDfRow")
    nCap   <- attr(x, "nCap")
    jMem   <- attr(x, "jMem")

    footerLines <- list()

    # Add truncated data message if applicable
    if (nDfRow > nCap) {
        footerLines[[length(footerLines) + 1]] <- pillar::style_subtle(sprintf(
            "# %s showing the first %s rows%s",
            cli::symbol$info,
            nCap,
            cli::symbol$ellipsis
        ))
    }

    # Always add Java memory address
    footerLines[[length(footerLines) + 1]] <- pillar::style_subtle(sprintf(
        "# %s Java memory address: 0x%s",
        cli::symbol$info,
        cli::style_bold(jMem)
    ))

    c(defaultFooter, unlist(footerLines))
}



# /// S4 - class definition /////////////////////////////////////////

## ----
#' @title
#' TasselPhenotype Class Definition
#'
#' @description
#' Defines the \code{TasselPhenotype} class, which represents
#' phenotype data in the TASSEL 5 framework.
#'
#' @slot attrData
#' A \code{tbl_df} containing attribute data for the phenotype.
#' @slot attrSummary
#' A \code{list} summarizing the attributes of the phenotype data.
#' @slot dispData
#' A \code{java_pheno_tbl} object for displaying phenotype data.
#' @slot rData
#' A \code{tbl_df} containing the phenotype data in R format.
#' @slot jRefObj
#' A \code{jobjRef} representing a reference to the Java object in
#' TASSEL 5.
#' @slot jMemAddress
#' A \code{character} string representing the memory address of the
#' Java object.
#' @slot jClass
#' A \code{character} string representing the Java class name of the
#' object.
#'
#' @name TasselPhenotype-class
#' @rdname TasselPhenotype-class
#' @exportClass TasselPhenotype
setClass(
    Class = "TasselPhenotype",
    slots = c(
        attrData    = "data.frame",
        attrSummary = "list",
        dispData    = "java_pheno_tbl",
        rData       = "data.frame",
        jRefObj     = "jobjRef",
        jMemAddress = "character",
        jClass      = "character"
    )
)


## ----
#' @title
#' Read and convert phenotype data into TASSEL 5 phenotype objects
#'
#' @description
#' This function reads phenotype data from either a file path or a
#' data frame.
#'
#' @details
#' \itemize{
#'   \item
#'   If \code{x} is a character string, the function assumes it is
#'   a file path and calls \code{readPhenotypeFromFile(x)}.
#'
#'   \item
#'   If \code{x} is a data frame, the function requires the
#'   \code{attr} parameter to provide metadata and calls
#'   \code{readPhenotypeFromDf(x, attr)}.
#'
#'   \item
#'   If \code{x} is neither a character string nor a data frame,
#'   the function throws an error.
#'
#' }
#'
#' @param x
#' A character string representing the file path to the phenotype data
#' or a data frame containing the phenotype data.
#' @param attr
#' An optional attribute metadata parameter required when \code{x} is
#' a data frame. Defaults to \code{NULL}.
#'
#' @return A phenotype object created from the input data.
#'
#' @examples
#' \dontrun{
#' # Reading phenotype data from a file
#' phenotype <- readPhenotype("path/to/phenotype/file.txt")
#'
#' # Reading phenotype data from a data frame
#' attrDf <- tibble::tribble(
#'     ~"col_id",      ~"tassel_attr",
#'     "taxa_id",      "taxa",
#'     "plant_height", "data",
#'     "PC1",          "covariate",
#'     "yield",        "data",
#' )
#' df <- tibble::tribble(
#'     ~"taxa_id", ~"plant_height", ~"PC1", ~"yield",
#'     "line_a",   12.3,            0.5,    2,
#'     "line_b",   22.8,            -1.5,   3,
#' )
#'
#' phenotypeDf <- readPhenotype(df, attr = attrDf)
#' }
#'
#' @export
readPhenotype <- function(x, attr = NULL) {
    if (is.character(x)) {
        return(readPhenotypeFromFile(x))
    } else if (is.data.frame(x)) {
        if (is.null(attr)) {
            rlang::abort("Phenotype objects evaluated from 'data.frame' need attribute metadata ('attr' parameter)")
        }
        return(readPhenotypeFromDf(x, attr))
    } else {
        rlang::abort("Unsupported input type for 'x'. Must be a file path ('character') or 'data.frame'")
    }
}



# /// Methods (show) ////////////////////////////////////////////////

## ----
#' @title
#' Display summary information of a TasselPhenotype object
#'
#' @param object
#' A \code{TasselPhenotype} object
setMethod("show", "TasselPhenotype", function(object) {
    print(object@dispData)
})



# /// Methods (general) /////////////////////////////////////////////

## ----
#' @rdname attributeData
#' @export
setMethod(
    f = "attributeData",
    signature = signature(object = "TasselPhenotype"),
    definition = function(object) {
        return(object@attrData)
    }
)


## ----
#' @rdname javaRefObj
#' @export
setMethod(
    f = "javaRefObj",
    signature = signature(object = "TasselPhenotype"),
    definition = function(object) {
        return(object@jRefObj)
    }
)


