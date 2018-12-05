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

## Source files
source("R/AllClasses.R")
source("R/AllGenerics.R")
source("R/TasselPluginWrappers.R")


## jinit
rJava::.jinit()

## Add class path
homeloc <- Sys.getenv("HOME")
rJava::.jaddClassPath(
    paste0(homeloc, "/Development/tassel_5_standalone/lib")
)
rJava::.jaddClassPath(
    paste0(homeloc, "/Development/tassel_5_standalone/sTASSEL.jar")
)


## Add VCF path
devloc <- paste0(homeloc, "/Development")
vcfPath <- paste0(
    devloc,
    "/tassel_5_test/dataFiles/GenotypeTableTests/",
    "maize_chr9_10thin100.recode.vcf"
)


