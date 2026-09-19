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
#'   \code{attrTypes} parameter to map each column to a TASSEL
#'   attribute type and calls
#'   \code{readPhenotypeFromDf(x, attrTypes)}.
#'
#'   \item
#'   If \code{x} is neither a character string nor a data frame,
#'   the function throws an error.
#'
#' }
#'
#' A file carries its own attribute types in its header, but a plain
#' data frame carries none, so \code{attrTypes} supplies them. It
#' plays the part \code{colData} plays in a
#' \code{SummarizedExperiment}: it must describe every column of
#' \code{x} exactly once, and describe nothing else. Nothing is read
#' off column position or guessed from the data, so a column left
#' out, named twice, or misspelled is an error rather than a silent
#' drop.
#'
#' Each attribute type accepts one R type, and the column is coerced
#' to what TASSEL needs to hold it:
#'
#' \describe{
#'   \item{\code{"taxa"}}{
#'     \code{character} or \code{factor}, and no missing values,
#'     since each observation has to name its taxon. Exactly one
#'     column is the taxa column.
#'   }
#'   \item{\code{"data"}, \code{"covariate"}}{
#'     \code{double} or \code{integer}. \code{NA} and non-finite
#'     values are carried through as TASSEL missing values.
#'   }
#'   \item{\code{"factor"}}{
#'     \code{character}, \code{factor}, \code{integer}, or
#'     \code{logical}, and no missing values, since TASSEL would read
#'     them as a category of their own. A \code{double} column is
#'     rejected, as rounding a continuous measure into categories is
#'     rarely intended.
#'   }
#' }
#'
#' Any other column type, such as a \code{Date} or a list column, is
#' rejected by name.
#'
#' @param x
#' A character string representing the file path to the phenotype data
#' or a data frame containing the phenotype data.
#' @param attrTypes
#' A mapping of each column of \code{x} to its TASSEL attribute type,
#' required when \code{x} is a data frame and ignored otherwise.
#' Either a named character vector, where the names are column names,
#' or a data frame. A data frame is read in either spelling:
#' \code{col_id} and \code{tassel_attr}, or the \code{trait_id} and
#' \code{trait_type} that \code{\link{attributeData}()} reports, so the
#' metadata read off one phenotype can be used to build another.
#' @param attr
#' \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
#' Renamed to \code{attrTypes}.
#'
#' @return A phenotype object created from the input data.
#'
#' @examples
#' \dontrun{
#' # Reading phenotype data from a file
#' phenotype <- readPhenotype("path/to/phenotype/file.txt")
#'
#' # Reading phenotype data from a data frame
#' df <- tibble::tribble(
#'     ~"taxa_id", ~"plant_height", ~"PC1", ~"yield",
#'     "line_a",   12.3,            0.5,    2,
#'     "line_b",   22.8,            -1.5,   3,
#' )
#'
#' phenotypeDf <- readPhenotype(
#'     df,
#'     attrTypes = c(
#'         taxa_id      = "taxa",
#'         plant_height = "data",
#'         PC1          = "covariate",
#'         yield        = "data"
#'     )
#' )
#'
#' # The same mapping, written as a data frame
#' attrDf <- tibble::tribble(
#'     ~"col_id",      ~"tassel_attr",
#'     "taxa_id",      "taxa",
#'     "plant_height", "data",
#'     "PC1",          "covariate",
#'     "yield",        "data",
#' )
#' readPhenotype(df, attrTypes = attrDf)
#'
#' # The return trip, using the metadata of an existing phenotype
#' readPhenotype(
#'     as.data.frame(phenotype),
#'     attrTypes = attributeData(phenotype)
#' )
#' }
#'
#' @export
readPhenotype <- function(x, attrTypes = NULL, attr = lifecycle::deprecated()) {
    if (lifecycle::is_present(attr)) {
        lifecycle::deprecate_warn(
            "0.14.0", "readPhenotype(attr)", "readPhenotype(attrTypes)"
        )
        if (is.null(attrTypes)) attrTypes <- attr
    }

    if (is.character(x)) {
        return(readPhenotypeFromFile(x))
    } else if (is.data.frame(x)) {
        if (is.null(attrTypes)) {
            rlang::abort(c(
                "A `data.frame` phenotype needs attribute metadata (`attrTypes`)",
                "i" = paste0(
                    "Pass a named character vector, e.g. ",
                    "c(Taxon = \"taxa\", EarHT = \"data\"), or a data frame with ",
                    "`col_id` and `tassel_attr` columns."
                )
            ))
        }
        return(readPhenotypeFromDf(x, attrTypes))
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



# /// Bracket Method /////////////////////////////////////////////////

## ----
#' @title Subset a TasselPhenotype
#'
#' @description
#' Matrix-style subsetting for \code{TasselPhenotype} objects using
#' the \code{ph[observations, traits]} syntax, the phenotype
#' counterpart of the \code{gt[taxa, sites]} syntax used on genotype
#' tables.
#'
#' @details
#' The row axis is observations rather than taxa, since a phenotype
#' may hold several observations of one taxon.
#' \code{\link{taxaWhere}()} tests each observation on its own, while
#' \code{\link{taxa}()} and a bare character vector keep whole taxa
#' with every observation they have.
#'
#' The column axis is traits. The taxa column is the label on the row
#' axis rather than a trait, so it is not a position in \code{j} and
#' is never dropped.
#'
#' Indices are applied left to right, so a trait predicate sees the
#' observations that the row index left behind. \code{notMissing} in
#' \code{\link{traitsWhere}()} is therefore computed over the
#' surviving observations, as it is when \code{\link{filterTraits}()}
#' follows \code{\link{filterTaxa}()} in a pipeline.
#'
#' @param x A \code{TasselPhenotype} object.
#' @param i Observation selector: a character vector of taxa IDs, a
#'   \code{\linkS4class{TaxaSelector}}, or missing.
#' @param j Trait selector: an integer vector of 1-based trait
#'   positions, a character vector of trait names, a
#'   \code{\linkS4class{TraitSelector}}, or missing.
#' @param ... Ignored.
#' @param drop Ignored.
#'
#' @return A new \code{TasselPhenotype} containing the selected
#'   observations and/or traits.
#'
#' @examples
#' \dontrun{
#' ph[taxa("33-16", "38-11"), ]
#' ph[taxaWhere(EarHT > 100), ]
#' ph[, traits("EarHT", "dpoll")]
#' ph[, traitsWhere(traitType == "covariate")]
#' ph[taxaWhere(location == "A"), 1:3]
#' ph[, !traits("EarDia")]
#' }
#'
#' @rdname TasselPhenotype-class
#' @aliases [,TasselPhenotype,ANY,ANY-method
setMethod("[", "TasselPhenotype", function(x, i, j, ..., drop = FALSE) {
    tasIn <- .resolveTasselInput(x, "phenotype", "[")

    jPh   <- x@jRefObj
    rData <- x@rData

    if (!missing(i)) {
        kept  <- applyPhenotypeTaxaSelector(jPh, i, tasIn, rData)
        jPh   <- kept$jPh
        rData <- rData[kept$positions, , drop = FALSE]
    }

    # An observation subset keeps every attribute, so the attribute
    # metadata cached on 'x' still describes the traits of 'jPh'
    if (!missing(j)) jPh <- applyTraitSelector(jPh, j, x@attrData, rData)

    createTasselPhenotype(jPh)
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


## ----
#' @rdname traitNames
#' @aliases traitNames,TasselPhenotype-method
#' @export
setMethod(
    f = "traitNames",
    signature = signature(object = "TasselPhenotype"),
    definition = function(object) {
        attrData <- object@attrData

        return(attrData$trait_id[attrData$trait_type != "taxa"])
    }
)



# /// Methods (taxa) ////////////////////////////////////////////////

## ----
#' @rdname taxaList
#' @aliases taxaList,TasselPhenotype-method
#' @export
setMethod("taxaList", "TasselPhenotype", function(tasObj) {
    .taxaNames(tasObj@jRefObj$taxa())
})



# /// Methods (coercion) ////////////////////////////////////////////

## ----
#' @title Coerce phenotype data to a data frame
#'
#' @description
#' Returns the phenotype data held by a \code{TasselPhenotype} object as a
#' \code{tibble}. Attribute metadata is available separately via
#' \code{\link{attributeData}()}.
#'
#' @param x A \code{TasselPhenotype} object.
#' @param row.names Ignored, present for generic compatibility.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Additional arguments to be passed to or from methods.
#'
#' @return A \code{tibble} of phenotype data.
#'
#' @export
as.data.frame.TasselPhenotype <- function(
    x,
    row.names = NULL,
    optional = FALSE,
    ...
) {
    return(x@rData)
}


