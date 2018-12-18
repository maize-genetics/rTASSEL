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

# Preamble

## Load packages
library(rJava) 
library(GenomicRanges)
library(stringr)
library(SummarizedExperiment)
library(snpStats)
library(hexbin)


## Set WD
setwd("~/Projects/rtassel/")

## Set Java jar path
path_tassel <- paste0(getwd(),"/inst/java/sTASSEL.jar")

## jinit
rJava::.jinit(parameters="-Xmx6g")
.jcall(.jnew("java/lang/Runtime"), "J", "totalMemory")
.jcall(.jnew("java/lang/Runtime"), "J", "maxMemory")

## Add class path
# Note the file class paths may differ between Windows and Macs.
homeloc <- Sys.getenv("HOME")

rJava::.jaddClassPath(path_tassel)
print(.jclassPath())

tasselVersion <- rJava::.jfield("net/maizegenetics/tassel/TASSELMainFrame","S","version")
str_c("Using TASSEL version: ",tasselVersion)

## Source files
source("R/AllGenerics.R")
source("R/AllClasses.R")
source("R/TasselPluginWrappers.R")
source("R/PullFunctions.R")
source("R/PushFunctions.R")
source("R/gwasPolyObjectCreator.R")

## Add VCF path
vcfPath <- paste0(
    getwd(),
    "/data/maize_chr9_10thin40000.recode.vcf"
)



# Tests

## Make genotype table
tasGenoTable <- readGenotypeTable(vcfPath)

## Make sample data frame of Taxa
tasDF <- sampleDataFrame(tasGenoTable)

## Make genomic ranges
tasGRanges <- genomicRanges(tasGenoTable)
tasGRanges
