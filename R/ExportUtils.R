#--------------------------------------------------------------------
# Script Name:   ExportUtils.R
# Description:   Export utilities for rTASSEL
# Author:        Brandon Monier
# Created:       2020-11-09 at 11:53:02
# Last Modified: 2020-11-10 at 17:02:40
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house methods for
#    exporting various rTASSEL data objects to flat files.
#--------------------------------------------------------------------

## ----
#' @title Export Genotype Table to Disk
#'
#' @description Exports genotype tables to various flat file formats.
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype} that
#'   contains a genotype table.
#' @param file Output file name.
#' @param format Export file format.
#' @param keepDepth Whether to keep depth if format supports depth. Defaults
#'   to \code{TRUE}.
#' @param taxaAnnotations Whether to include taxa annotations if format
#'   supports taxa. Defaults to \code{TRUE}.
#' @param branchLengths Whether to include branch lengths for Newick formatted
#'   files. Defaults to \code{TRUE}.
#'
#' @importFrom rJava .jchar
#' @importFrom rJava is.jnull
#' @importFrom rJava J
#'
#' @export
exportGenotypeTable <- function(tasObj,
                                file = "",
                                format = c("vcf", "hapmap", "plink", "flapjack", "hdf5"),
                                keepDepth = TRUE,
                                taxaAnnotations = TRUE,
                                branchLengths = TRUE) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    # Filter type selection
    acceptedFormats <- c("vcf", "hapmap", "plink", "flapjack", "hdf5")
    format <- match.arg(format)
    if (missing(format) || !format %in% acceptedFormats) {
        stop(
            paste(
                "Please specify analysis type",
                "(\"vcf\", \"hapmap\", \"plink\", \"HDF5\", or \"flapjack\")"
            )
        )
    }

    rJC <- rJava::J("net.maizegenetics.dna.snp.ExportUtils")

    if (format == "vcf") {
        rJC$writeToVCF(
            jGenoTable,
            file,
            keepDepth, # <- keep depth
            NULL  # <- progress listener
        )
    } else if (format == "hapmap") {
        rJC$writeToHapmap(
            jGenoTable,
            FALSE, # <- diploid?
            file,
            rJava::.jchar(9), # <- tab delimited UTF char
            taxaAnnotations,
            NULL
        )
    } else if (format == "plink") {
        rJC$writeToPlink(
            jGenoTable,
            file,
            rJava::.jchar(9) # <- tab delimited UTF char
        )
    } else if (format == "hdf5") {
        rJC$writeGenotypeHDF5(
            jGenoTable,
            file,
            NULL, # <- export taxa annotation
            keepDepth
        )
    } else {
        message("FlapJack format not yet implemented. Coming soon!")
    }
}


