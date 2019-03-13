#--------------------------------------------------------------------
# Script Name:   GLMFunctions.R
# Description:   GLM related functions
# Author:        Brandon Monier
# Created:       2019-03-06 at 17:50:59
# Last Modified: 2019-03-12 at 17:23:49
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house various functions
#    necessary for GLM operations for TASSEL
#--------------------------------------------------------------------

# Run TASSEL GLM on Genotype Phenotype object
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
    parser <- function(obj) {
        obj <- unlist(strsplit(obj$toStringTabDelim(), split = "\n"))
        obj <- strsplit(obj, split = "\t")
        obj <- t(simplify2array(obj))
        colnames(obj) <- as.character(unlist(obj[1, ]))
        obj <- obj[-1, ]
        tibble::as_tibble(obj)
    }

    glm_stats <- parser(glm_stats)
    glm_geno  <- parser(glm_geno)

    glmColConvert(
        stat = glm_stats,
        geno = glm_geno
    )

}

# Convert GLM columns to proper R data types
glmColConvert <- function(stat, geno) {
    # Numeric convert
    stat[4:18] <- sapply(stat[4:18], as.numeric)
    geno[c(4, 5, 7)] <- lapply(geno[c(4, 5, 7)], as.numeric)

    # Factor convert
    stat[c(1, 3)] <- lapply(stat[c(1, 3)], factor)
    geno[c(1, 3, 6)] <- lapply(geno[c(1, 3, 6)], factor)

    # Reorder Chromsome
    stat$Chr <- factor(
        stat$Chr,
        levels = paste(sort(as.numeric(levels(stat$Chr))))
    )
    geno$Chr <- factor(
        geno$Chr,
        levels = paste(sort(as.numeric(levels(geno$Chr))))
    )

    return(
        list(
            "GLM_Statistics" = tibble::as_tibble(stat),
            "GLM_Geno_Effects" = tibble::as_tibble(geno)
        )
    )
}
