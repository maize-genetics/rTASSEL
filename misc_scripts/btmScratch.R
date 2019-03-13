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
