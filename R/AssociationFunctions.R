#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   AssociationFunctions.R
# Description:   General functions for running association analysis
# Author:        Brandon Monier
# Created:       2019-03-29 at 13:51:02
# Last Modified:
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house all necessary
#    general functions and wrappers for TASSEL association analyses.
#--------------------------------------------------------------------

# Association analysis front-end
assocModelFitter <- function(tasObj, model, kinship = NULL) {

}



# TASSEL Table to R data frame converter
tasTableConvert <- function(stringTab) {
    obj <- unlist(strsplit(stringTab, split = "\n"))
    obj <- strsplit(obj, split = "\t")
    obj <- t(simplify2array(obj))
    colnames(obj) <- as.character(unlist(obj[1, ]))
    obj <- obj[-1, ]
    tibble::as_tibble(obj)
}
