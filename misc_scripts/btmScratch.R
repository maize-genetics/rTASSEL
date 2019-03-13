#                     _       _                     _
#                    | |     | |                   | |
#  ___  ___ _ __ __ _| |_ ___| |__  _ __   __ _  __| |
# / __|/ __| '__/ _` | __/ __| '_ \| '_ \ / _` |/ _` |
# \__ \ (__| | | (_| | || (__| | | | |_) | (_| | (_| |
# |___/\___|_|  \__,_|\__\___|_| |_| .__/ \__,_|\__,_|
#                                  | |
#                                  |_|

#=== Notes ==========================================================



#=== Functions ======================================================

## Association output for Terry
terryAssoc <- function(assocDF) {
    traitNames <- as.vector(assocDF[[1]])
    jc <- J("net/maizegenetics/plugindef/GenerateRCode")
    jc$association(traitNames)
}



# === Miscellaneous =================================================

## Create a manhattan plot
glmStats <- rTASSEL:::tasselGLM(tasGenoPhenoObj = tasGenoPheno)
glmStats <- glmStats$GLM_Statistics

ggplot2::ggplot(data = glmStats[glmStats$Trait == "dpoll", ]) +
    ggplot2::aes(x = Pos, y = -log10(p)) +
    ggplot2::facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(
            angle = 90, hjust = 1
        ),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank()
    ) +
    ggplot2::xlab("Position") +
    ggplot2::ylab(bquote(~-log[10]~ '('*italic(p)*'-value)'))






## `assocModelDesign()` debug - DON'T RUN
jtsPheno <- rTASSEL:::getPhenotypeTable(tasGenoPheno)
phenoAttDf <- rTASSEL:::extractPhenotypeAttDf(jtsPheno)
phenoAttDf <- tibble::add_case(
    phenoAttDf,
    traitName = "G",
    traitType = "genotype",
    traitAttribute = "Genotype"
)
df <- rTASSEL:::emptyDFWithPhenotype(phenoAttDf)

