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
#' @name TasselPhenotype-class
#' @rdname TasselPhenotype-class
#' @exportClass TasselPhenotype
setClass(
    Class = "TasselPhenotype",
    slots = c(
        attrData    = "tbl_df",
        attrSummary = "list",
        dispData    = "java_pheno_tbl",
        rData       = "tbl_df",
        jPheno      = "jobjRef",
        jMemAddress = "character",
        jClass      = "character"
    )
)


## ----
#' @title Helper function to build a TasselPhenotype object
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

#' @export
setMethod("show", "TasselPhenotype", function(object) {
    print(object@dispData)
})



# /// Methods (general) /////////////////////////////////////////////


