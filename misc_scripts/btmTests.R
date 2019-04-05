#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   btmTests.R
# Description:   GenomicRanges test for rTASSEL
# Author:        Brandon Monier
# Created:       2018-12-04 at 16:46:06
# Last Modified: 2019-03-13 at 15:07:04
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to test out wrapper
#    capabilities for rTASSEL in generating GenomicRanges class
#--------------------------------------------------------------------

# === Preamble (just load rTASSEL instead) ==========================

# ## jinit memory test
# rJava::.jinit(parameters="-Xmx6g")
# .jcall(.jnew("java/lang/Runtime"), "J", "totalMemory")
# .jcall(.jnew("java/lang/Runtime"), "J", "maxMemory")

## Start logging file
rTASSEL::startLogger(fullPath = "/home/bm646/Temporary/")



# === Get paths =====================================================

## Genotype file path example
genoPath <- system.file("extdata", "mdp_genotype.hmp.txt", package = "rTASSEL")

## Phenotype file path example
phenoPath  <- system.file("extdata", "mdp_traits.txt", package = "rTASSEL")
phenoPath2 <- system.file("extdata", "mdp_phenotype.txt", package = "rTASSEL")
phenoPath3 <- system.file("extdata", "mdp_traits_nomissing.txt", package = "rTASSEL")



# === Get phenotype data frame ======================================

## Phenotype data frame example - currently not working
phenoDF <- read.table(phenoPath, header = TRUE)
colnames(phenoDF)[1] <- "Taxon"



# === Read tests ====================================================

## Read GenotypeTable
tasGeno <- rTASSEL::readGenotypeTableFromPath(genoPath)

## Read PhenotypeTable
tasPheno <- rTASSEL::readPhenotypeFromPath(phenoPath)
tasPheno2 <- rTASSEL::readPhenotypeFromPath(phenoPath2)
tasPheno3 <- rTASSEL::readPhenotypeFromDataFrame(
    phenotypeDF = phenoDF,
    taxaID = "Taxon",
    attributeTypes = c("data", "covariate", "data")
)

## Read Genotype and Phenotype
tasGenoPheno <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = genoPath,
    phenoPathDFOrObj = phenoPath
)

## Read Genotype and Phenotype (with covariates)
tasGenoPhenoCov <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = genoPath,
    phenoPathDFOrObj = phenoPath2
)

## Read Genotype and Phenotype (no missing values)
tasGenoPhenoFast <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = genoPath,
    phenoPathDFOrObj = phenoPath3
)

## Read Genotype and Phenotype (from two `TasselGenotypePhenotype` objects)
tasGenoPhenoTASOBJ <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = tasGeno,
    phenoPathDFOrObj = tasPheno
)

## Read Genotype and Phenotype (from Genotype path and R phenotype data frame)
tasGenoPhenoDF <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = genoPath,
    phenoPathDFOrObj = phenoDF,
    taxaID = "Taxon"
)


# === Test SummarizedExperiment from TASSEL Genotype ================
tasSumExp <- rTASSEL::getSumExpFromGenotypeTable(tasGenoPheno)
tasGRanges <- GenomicRanges::granges(tasSumExp)



# === Test extract phenotype data frame =============================
tasTmpPheno <- rTASSEL::getPhenotypeDF(tasGenoPheno)
tasTmpPhenoCov <- rTASSEL::getPhenotypeDF(tasGenoPhenoCov)



# === Matrix functions =============================================+

## Test `kinshipMatrix()` - return kinship matrix TASSEL object
##     calculated from a `TasselGenotypePhenotype` object
tasKin <- rTASSEL::kinshipMatrix(
    tasObj = tasGenoPheno,
    method = "Centered_IBS",
    maxAlleles = 6,
    algorithmVariation = "Observed_Allele_Freq"
)
tasKinR <- rTASSEL::kinshipToRMatrix(tasKin) ## Convert to R matrix
tasKinR[1:10, 1:10]                        ## Get subset
image(tasKinR)                             ## Visualize it


## Test `distanceMatrix()` - return distance matrix TASSEL object
tasDist <- rTASSEL::distanceMatrix(tasGenoPheno)
tasDistR <- rTASSEL::distanceToRMatrix(tasDist)
tasDistR[1:10, 1:10]
image(tasDistR)


# === Test association functions ====================================

tasBLUE <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPhenoCov,
    formula = list(EarHT, dpoll) ~ .,
    fitMarkers = FALSE,
    kinship = NULL,
    fastAssociation = FALSE
); tasBLUE

tasGLM <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPhenoCov,
    formula = list(EarHT, dpoll) ~ .,
    fitMarkers = TRUE,
    kinship = NULL,
    fastAssociation = FALSE
); tasGLM

tasMLM <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPheno,
    formula = . ~ .,
    fitMarkers = TRUE,
    kinship = tasKin,
    fastAssociation = FALSE
); tasMLM

tasFAST <- rTASSEL::assocModelFitter(
    tasObj = tasGenoPhenoFast,
    formula = . ~ .,
    fitMarkers = TRUE,
    kinship = NULL,
    fastAssociation = TRUE
); tasFAST


## Test WRONG file path (readGenotypeTable())
# genoWrongPath <- "this/is/so/wrong"
# rTASSEL::readGenotypeTable(path = genoWrongPath)

## Test WRONG file path (readPhenotypeTable())
# phenoWrongPath <- "wrong/path"
# rTASSEL::readPhenotypeTable(path = phenoWrongPath)

## Test GLM functionality
# tasGLM <- rTASSEL:::tasselGLM(tasGenoPheno)



# === Visualizations for association data ===========================

## Test Manhattan plot functionality
rTASSEL:::manhattanPlot(
    assocStats = tasGLM$GLM_Statistics,
    trait = "dpoll",
    threshold = 10
)
