#--------------------------------------------------------------------
# Script Name:   RAmpSeq.R
# Description:   rAmpSeq pipeline
# Author:        Brandon Monier
# Created:       2021-02-17 at 21:50:29
# Last Modified: 2021-02-17 at 22:14:11
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house methods for
#    running the rAmpSeq pipeline using rTASSEL instead of
#    compiling the Java source code with static fields.
#--------------------------------------------------------------------


#' @title Run rAmpSeq in R
#'
#' @description Runs the rAmpSeq pipeline in R with user defined methods.
#'
#' @param refGenomeFile A reference genome file (e.g. \code{*.fa.gz}).
#' @param kmerLength Length of k-mer to be identified.
#' @param minAmpLength Minimum length of amplicons to be produced.
#' @param maxAmpLength Maximum length of amplicons to be produced.
#' @param maxPrimerObs Maximum number of "good" primers to be selected.
#' @param minCopyNum Minimum amplicon copy threshold.
#' @param maxCopyNum Maximum amplicon copy threshold.
#' @param ambigBPs Do you want ambiguous basepairs to be generated? Defaults
#'   to \code{FALSE}. See details for further explanation.
#'
#' @details Currently, the \code{ambigBPs} parameters convert every third
#'   base pair to an \code{A}. For example, if a k-mer of length 12 was
#'   generated:
#'
#'   \code{A C C G T A T T C G A C}
#'
#'   would be converted to:
#'
#'   \code{A C \strong{A} G T \strong{A} T T \strong{A} G A \strong{A}}
#'
#' @importFrom rJava .jnew
#'
#' @export
rAmpSeq <- function(
    refGenomFile,
    kmerLength,
    minAmpLength,
    maxAmpLenght,
    maxPrimerObs,
    minCopyNum,
    maxCopyNum,
    ambigBPs = FALSE
) {
    rJava::.jnew(
        class = "com.btmonier.rampseq.RAmpSeq",
        refGenomeFile,
        as.integer(kmerLength),
        as.integer(minAmpLength),
        as.integer(maxAmpLength),
        as.integer(maxPrimerObs),
        as.integer(minCopyNum),
        as.integer(maxCopyNum),
        ambigBPs,
        FALSE,
        FALSE, FALSE, FALSE, FALSE, # positional degeneracy
        "kmerMap.ser",
        "kmerPosMap.ser"
    )
}


