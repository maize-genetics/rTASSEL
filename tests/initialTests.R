#Test script for R classes and functions

##Installs needed
# install.packages("tidyverse")
# To install Java on a mac
# make sure latest git/gcc tools are installed
# from terminal: xcode-select --install
# With R closed, run this to let R know what Java is available.  Close the term when done
# from terminal: R CMD javareconf
# install.packages("rJava", dependencies=TRUE, type="source")
# if (!requireNamespace("BiocManager")) {
#   install.packages("BiocManager")
# }
# BiocManager::install("SummarizedExperiment")

## Load packages
library(rJava) 
library(GenomicRanges)
library(stringr)


is_experimental <- TRUE
## Set WD
setwd("~/Code/rtassel")
path_exp_tassel <- paste0(getwd(),"/inst/java/sTASSEL.jar")



## jinit
rJava::.jinit()
.jcall(.jnew("java/lang/Runtime"), "J", "totalMemory")
.jcall(.jnew("java/lang/Runtime"), "J", "maxMemory")

## Add class path
# Note the file class paths may differ between Windows and Macs.
homeloc <- Sys.getenv("HOME")
if(is_experimental == TRUE) {
  tasselPath <- path_exp_tassel
} else {
  tasselPath <- paste0(homeloc, "/Code/tassel-5-standalone/sTASSEL.jar")
}

rJava::.jaddClassPath(tasselPath)
print(.jclassPath())

tasselVersion <- rJava::.jfield("net/maizegenetics/tassel/TASSELMainFrame","S","version")
str_c("Using TASSEL version: ",tasselVersion)

# rJava::.jaddClassPath("/Users/edwardbuckler/Code/tassel-5-source/dist/sTASSEL.jar")

## Source files
source("R/AllGenerics.R")
source("R/AllClasses.R")
source("R/TasselPluginWrappers.R")
source("R/PullFunctions.R")


## VCF file path exmaple...
vcfPath <- paste0(
  getwd(),
  "/data/maize_chr9_10thin40000.recode.vcf"
)

test <- readGenotypeTable(vcfPath)
test
test@name
test@jtsGenotypeTable
test2 <- filterSiteBuilderPlugin(test, siteMinCount = 0)
test2
test2@name
test2@jtsGenotypeTable
testPostions <- positions(test)
testPostions
testTaxa <- taxa(test)
testTaxa

sampleDF <- sampleDataFrame(testTaxa)

genoIntArray <- rJava::.jcall(
  "net/maizegenetics/plugindef/GenerateRCode",
  "[I",
  "genotypeTableToDosageIntArray",
  test@jtsGenotypeTable
)

genoNameArray <- rJava::.jcall(
  "net/maizegenetics/plugindef/GenerateRCode",
  "[S",
  "genotypeTableToSampleNameArray",
  test@jtsGenotypeTable
)

genoPositionVector <- rJava::.jcall(
  "net/maizegenetics/plugindef/GenerateRCode",
  "Lnet/maizegenetics/plugindef/GenerateRCode$PositionVectors;",
  "genotypeTableToPositionListOfArrays",
  test@jtsGenotypeTable
)

