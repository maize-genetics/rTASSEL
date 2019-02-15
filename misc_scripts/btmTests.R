#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   btmTests.R
# Description:   GenomicRanges test for rTASSEL
# Author:        Brandon Monier
# Created:       2018-12-04 at 16:46:06
# Last Modified: 2019-01-23 at 14:36:19
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to test out wrapper
#    capabilities for rTASSEL in generating GenomicRanges class
#--------------------------------------------------------------------

# === Create log file and output messages from console ==============

## NOTE: For my use only (START AT PREAMBLE)
setwd("~/Projects/rtassel")
path_tassel <- paste0(getwd(),"/inst/java/sTASSEL.jar")
rJava::.jinit(parameters="-Xmx6g")
rJava::.jaddClassPath(path_tassel)

## Make rtassel_output file
rtOut <- paste0(Sys.getenv("HOME"), "/Temporary/rtassel_output")
if (!exists(rtOut)) {
    system(paste("touch", rtOut))
}

## Send TASSEL console output messages to file
rJava::.jcall(
    "net.maizegenetics/util/LoggingUtils",
    "V",
    "setupLoggingOff"
)



# === Preamble ======================================================

## Set WD (local use only)
setwd("~/Projects/rtassel")

## Load packages
library(rJava)
library(GenomicRanges)
library(stringr)
library(SummarizedExperiment)
library(snpStats)
library(hexbin)
library(S4Vectors)

## jinit
rJava::.jinit(parameters="-Xmx6g")
.jcall(.jnew("java/lang/Runtime"), "J", "totalMemory")
.jcall(.jnew("java/lang/Runtime"), "J", "maxMemory")

## Add class path
path_tassel <- "inst/java/sTASSEL.jar"
rJava::.jaddClassPath(path_tassel)

## Which Tassel version
tasselVersion <- rJava::.jfield("net/maizegenetics/tassel/TASSELMainFrame","S","version")
paste0("Using TASSEL version: ", tasselVersion)

## Source files
source("R/AllGenerics.R")
source("R/AllClasses.R")
source("R/GenotypePhenotypeFunctions.R")
source("R/TasselPluginWrappers.R")
source("R/GenotypeTableFunctions.R")
source("R/PhenotypeFunctions.R")
source("R/PositionListFunctions.R")
source("R/TaxaListFunctions.R")
source("R/TableReportFunctions.R")
source("R/PluginSupport.R")



# === Tests =========================================================

## Genotype file path example
genoPath <- "inst/extdata/mdp_genotype.hmp.txt"

## Phenotype file path example
phenoPath <- "inst/extdata/mdp_traits.txt"
phenoPath2 <- "inst/extdata/mdp_phenotype.txt"

## Phenotype data frame example
phenoDF <- read.table(phenoPath, header = TRUE)
colnames(phenoDF)[1] <- "Taxon"

## Read GenotypeTable
tasGeno <- readGenotypeTable(genoPath)

## Read PhenotypeTable
tasPheno <- readPhenotypeTable(phenoPath)
tasPheno2 <- readPhenotypeTable(phenoPath2)

## Read Genotype and Phenotype
tasGenoPheno <- readGenotypePhenotype(
    genoPathOrObj = genoPath,
    phenoPathDFOrObj = phenoPath2
)

## Test `assocModelDesign()` - should return BLUE option
assocModelDesign(
    phenotypeGenotype = tasGenoPheno,
    fixed = list(EarHT, EarDia) ~ location + Q1 + Q2 + Q3 + Taxa,
    kinship = NULL
)

## Test `assocModelDesign()` - should return GLM option
assocModelDesign(
    phenotypeGenotype = tasGenoPheno,
    fixed = list(EarHT, EarDia) ~ location + Q1 + Q2 + Q3 + G,
    kinship = NULL
)

## Test `assocModelDesign()` - should return MLM option
assocModelDesign(
    phenotypeGenotype = tasGenoPheno,
    fixed = list(EarHT, dpoll) ~ location + G,
    kinship = "K"
)

## Test `kinshipPlugin()` - return kinship matrix TASSEL object
##     calculated from a `TasselGenotypePhenotype` object
tasKin <- kinshipPlugin(
    genotypeTable = tasGenoPheno,
    method = "Centered_IBS",
    maxAlleles = 6,
    algorithmVariation = "Observed_Allele_Freq"
)
tasKinR <- kinshipRMatrix(tasKin) ## Convert to R matrix
tasKinR[1:10, 1:10]               ## Get subset
image(tasKinR)                    ## Visualize it



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
