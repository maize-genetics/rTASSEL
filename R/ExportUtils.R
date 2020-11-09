#--------------------------------------------------------------------
# Script Name:   ExportUtils.R
# Description:   Export utilities for rTASSEL
# Author:        Brandon Monier
# Created:       2020-11-09 at 11:53:02
# Last Modified: 2020-11-09 at 12:12:29
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
#' @export
#'
exportGenotypeTable <- function(tasObj,
                                file = "",
                                format = c("vcf", "hapmap", "plink", "flapjack")) {
    # TODO get TASSEL methods
}
