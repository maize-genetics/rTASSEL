#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   btmTests.R
# Description:   GenomicRanges test for rTASSEL
# Author:        Brandon Monier
# Created:       2018-12-04 at 16:46:06
# Last Modified: 2018-12-04 at 19:00:38
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to test out wrapper
#    capabilities for rTASSEL in generating GenomicRanges class
#--------------------------------------------------------------------

# Create log file and output messages from console
if (!exists("~/Temporary/rtassel_output")) system("touch ~/Temporary/rtassel_output")
rJava::.jcall("net.maizegenetics/util/LoggingUtils", "V", "setupLogfile", "/home/bm646/Temporary/rtassel_output")


# Preamble

## Load packages
library(rJava) 
library(GenomicRanges)
library(stringr)
library(SummarizedExperiment)
library(snpStats)
library(hexbin)
library(S4Vectors)

## Set WD
setwd("~/Projects/rtassel")

## jinit
rJava::.jinit(parameters="-Xmx6g")
.jcall(.jnew("java/lang/Runtime"), "J", "totalMemory")
.jcall(.jnew("java/lang/Runtime"), "J", "maxMemory")

## Add class path
# Note the file class paths may differ between Windows and Macs.
homeloc <- Sys.getenv("HOME")
path_tassel <- paste0(getwd(),"/inst/java/sTASSEL.jar")
rJava::.jaddClassPath(path_tassel)
print(.jclassPath())

tasselVersion <- rJava::.jfield("net/maizegenetics/tassel/TASSELMainFrame","S","version")
str_c("Using TASSEL version: ",tasselVersion)

# rJava::.jaddClassPath("/Users/edwardbuckler/Code/tassel-5-source/dist/sTASSEL.jar")

## Source files
source("R/AllGenerics.R")
source("R/AllClasses.R")
source("R/TasselPluginWrappers.R")
source("R/PullFunctions.R")
source("R/PushFunctions.R")



# Tests

## Genotype file path example
genoPath <- paste0(
  getwd(),
  "/data/maize_chr9_10thin40000.recode.vcf"
)

## Phenotype file path example
phenoPath <- paste0(
    getwd(),
    "/data/mdp_traits.txt"
)

## Phenotype data frame example
phenoDF <- read.table(phenoPath, header = TRUE)
colnames(phenoDF)[1] <- "Taxon"

## Read GenotypeTable
tasGeno <- readGenotypeTable(genoPath)

## Read PhenotypeTable
tasPheno <- readPhenotypeTable(phenoPath)
