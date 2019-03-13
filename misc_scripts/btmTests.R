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

## Phenotype data frame example - currently not working
phenoDF <- read.table(phenoPath, header = TRUE)
colnames(phenoDF)[1] <- "Taxon"



# === Read tests ====================================================

## Read GenotypeTable
tasGeno <- rTASSEL::readGenotypeTable(genoPath)

## Read PhenotypeTable
tasPheno <- rTASSEL::readPhenotypeTable(phenoPath)
tasPheno2 <- rTASSEL::readPhenotypeTable(phenoPath2)

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



# === Assocatiation tests ===========================================

## Test `assocModelDesign()` - should return BLUE option
assocDF <- rTASSEL::assocModelDesign(
    phenotypeGenotype = tasGenoPhenoCov,
    fixed = list(EarHT, EarDia) ~ location + Q1 + Q2 + Q3 + Taxa,
    kinship = NULL
)

## Test `assocModelDesign()` - should return GLM option
## G = genotype
rTASSEL::assocModelDesign(
    phenotypeGenotype = tasGenoPhenoCov,
    fixed = list(EarHT, EarDia) ~ location + Q1 + Q2 + Q3 + G,
    kinship = NULL
)

## Test `assocModelDesign()` - should return MLM option
rTASSEL::assocModelDesign(
    phenotypeGenotype = tasGenoPhenoCov,
    fixed = list(EarHT, dpoll) ~ location + G,
    kinship = "K"
)



# ==== Kinship functions ============================================

## Test `kinshipPlugin()` - return kinship matrix TASSEL object
##     calculated from a `TasselGenotypePhenotype` object
tasKin <- rTASSEL::kinshipPlugin(
    genotypeTable = tasGenoPheno,
    method = "Centered_IBS",
    maxAlleles = 6,
    algorithmVariation = "Observed_Allele_Freq"
)
tasKinR <- rTASSEL::kinshipRMatrix(tasKin) ## Convert to R matrix
tasKinR[1:10, 1:10]                        ## Get subset
image(tasKinR)                             ## Visualize it



# === Test for errors ===============================================

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
