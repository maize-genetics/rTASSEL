#Test script for R classes and functions

setwd("/Users/edwardbuckler/Code/rtassel/")
setwd()

#compiles the classes needed
source("R/AllClasses.R")

## VCF file path exmaple...
vcfPath <- paste0(
  "/Users/edwardbuckler/Code/tassel-5-test/dataFiles/GenotypeTableTests/",
  "maize_chr9_10thin100.recode.vcf"
)

## .jcalls
tsGenotypeTable <- rJava::.jcall(
  "net/maizegenetics/dna/snp/ImportUtils",
  "Lnet/maizegenetics/dna/snp/GenotypeTable;",
  "read",
  vcfPath
)

test <- readGenotypeTable(vcfPath)
test
test@name
test@jtsGenotypeTable
test2 <- filterSiteBuilderPlugin(test@jtsGenotypeTable, siteMinCount = 40)
test2
test2@name
test2@jtsGenotypeTable
