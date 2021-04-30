#---------------------------------------------------------------------
# Description:   Support working with TASSEL FactorTables
# Author:        Brandon Monier & Terry Casstevens
# Created:       2021-04-30 at 10:07:53
# Last Modified: 2021-04-30 at 10:07:53
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions
#    necessary for reading in haplotype/factor datasets into R and
#    extracting data from TASSEL 6 objects.
#--------------------------------------------------------------------

#' @title Wrapper function of FactorTable class for haplotype data data
#'
#' @description This function is a wrapper for the
#'    \code{FactorTable} class. It is used for storing genotype
#'    information into a class object.
#'
#' @name readFactorTable
#' @rdname readFactorTable
#'
#' @param file A haplotype or factor data path (e.g. \code{*.vcf}).
#' @param verbose Show messages to console? Defaults to \code{TRUE}.
#'
#' @return Returns an object of \code{FactorTable} class.
#'
#' @importFrom rJava .jnew
#' @importFrom methods new
#' @export
readFactorTable <- function(file, verbose = TRUE) {

    if (!file.exists(file)) {
        stop("Cannot open file ", file, ": No such file or directory")
    }

    message("Reading VCF file...")
    factorBuilder <- rJava::.jnew("net/maizegenetics/dna/factor/io/BuilderFromHaplotypeVCF")

    methods::new(
        Class = "FactorTable",
        jFactorTable = factorBuilder$read(file)
    )
}


