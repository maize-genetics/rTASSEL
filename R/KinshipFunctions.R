#--------------------------------------------------------------------
# Script Name:   KinshipFunctions.R
# Description:   Kinship-related functions
# Author:        Brandon Monier
# Created:       2019-02-14 at 15:16:20
# Last Modified: 2019-02-14 at 15:20:58
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house various functions
#    necessary for Kinship matrix operations for TASSEL.
#--------------------------------------------------------------------

#' @title Convert TASSEL kinship object to an R matrix class
#'
#' @param kinJobj a TASSEL kinship object
#'
#' @export
kinshipRMatrix <- function(kinJobj) {
    tmp1 <- unlist(strsplit(kinJobj$toStringTabDelim(), split = "\n"))
    tmp2 <- strsplit(tmp1, split = "\t")
    tmp3 <- t(simplify2array(tmp2))
    colnames(tmp3) <- as.character(unlist(tmp3[1, ]))
    tmp3 <- tmp3[-1, ]
    matRow <- tmp3[, 1]
    tmp3 <- tmp3[, -1]
    tmp3 <- apply(tmp3, 2, as.numeric)
    rownames(tmp3) <- matRow
    return(tmp3)
}
