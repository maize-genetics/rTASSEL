#--------------------------------------------------------------------
# Script Name:   ExportUtils.R
# Description:   Export utilities for rTASSEL
# Author:        Brandon Monier
# Created:       2020-11-09 at 11:53:02
# Last Modified: 2020-11-09 at 12:27:54
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
#'
#' @importFrom rJava is.jnull
#'
#' @export
exportGenotypeTable <- function(tasObj,
                                file = "",
                                format = c("vcf", "hapmap", "plink", "flapjack")) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    # Filter type selection
    format <- match.arg(format)
    if (missing(format) || !format %in% c("vcf", "hapmap", "plink", "flapjack")) {
        stop(
            paste(
                "Please specify analysis type",
                "(\"vcf\", \"hapmap\", \"plink\", or \"flapjack\")"
            )
        )
    }


    rJC <- rJava::J("net.maizegenetics.dna.snp.ExportUtils")

    if (format == "vcf") {
        rJC$writeToVCF(
            tasObj,
            file,
            TRUE, # <- keep depth
            NULL  # <- progress listener
        )
    } else {
        message("Format not yet implemented. Coming soon!")
    }
}


