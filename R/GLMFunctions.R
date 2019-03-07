#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   GLMFunctions.R
# Description:   GLM related functions
# Author:        Brandon Monier
# Created:       2019-03-06 at 17:50:59
# Last Modified:
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house various functions
#    necessary for GLM operations for TASSEL
#--------------------------------------------------------------------

tasselGLM <- function(tasGenoPhenoObj) {
    plugin <- new(
        J("net.maizegenetics.analysis.association.FixedEffectLMPlugin"),
        .jnull(),
        FALSE
    )
    jtsDataSet <- createTasselDataSet(c(tasGenoPhenoObj@jTasselObj))
    glm <- plugin$processData(jtsDataSet)

    ## TODO (Brandon) - Iterate this process
    glm_stats <- glm$getData(0L)$getData()
    glm_geno  <- glm$getData(1L)$getData()

    # Parse GLM statistics
    glm_stats <- unlist(strsplit(glm_stats$toStringTabDelim(), split = "\n"))
    glm_stats <- strsplit(glm_stats, split = "\t")
    glm_stats <- t(simplify2array(glm_stats))
    colnames(glm_stats) <- as.character(unlist(glm_stats[1, ]))
    glm_stats <- glm_stats[-1, ]

    # Parse GLM Genotype Effects
    glm_geno <- unlist(strsplit(glm_geno$toStringTabDelim(), split = "\n"))
    glm_geno <- strsplit(glm_geno, split = "\t")
    glm_geno <- t(simplify2array(glm_geno))
    colnames(glm_geno) <- as.character(unlist(glm_geno[1, ]))
    glm_geno <- glm_geno[-1, ]

    return(
        list(
            "GLM_Statistics" = tibble::as_tibble(glm_stats),
            "GLM_Geno_Effects" = tibble::as_tibble(glm_geno)
        )
    )
}
