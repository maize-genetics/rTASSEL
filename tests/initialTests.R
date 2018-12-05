#Test script for R classes and functions

setwd("/Users/edwardbuckler/Code/rtassel/")
setwd()

#compiles the classes needed
source("R/AllGenerics.R")
source("R/AllClasses.R")
source("R/TasselPluginWrappers.R")

## VCF file path exmaple...
vcfPath <- paste0(
  "/Users/edwardbuckler/Code/tassel-5-test/dataFiles/GenotypeTableTests/",
  "maize_chr9_10thin100.recode.vcf"
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
