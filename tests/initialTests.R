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
# BiocManager::install("snpStats")
# BiocManager::install("hexbin")

## Load packages
library(rJava) 
library(GenomicRanges)
library(stringr)
library(SummarizedExperiment)
library(snpStats)
library(hexbin)
library(testthat)

## Set WD
setwd("~/Code/rtassel")
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

#Load a VCF file from disk and a R wrapped TASSEL GenotypeTable
aGenoTable <- readGenotypeTable(vcfPath)
show(aGenoTable)
expect_equal(aGenoTable@jtsGenotypeTable$numberOfSites(), 493)
expect_equal(aGenoTable@jtsGenotypeTable$numberOfTaxa(), 189)

#Filter a genotype table based minimum count
siteFiltGenoTable <- filterSiteBuilderPlugin(aGenoTable, siteMinCount = 40)
show(siteFiltGenoTable)
expect_equal(siteFiltGenoTable@jtsGenotypeTable$numberOfSites(), 389)
expect_equal(siteFiltGenoTable@jtsGenotypeTable$numberOfTaxa(), 189)


#Extracting the positions wrapper from a GenotypeTable
aPositions <- positions(aGenoTable)
aPositions
#Extract genomicRanges from GenotypeTable, then convert to data.frame
gr <- genomicRanges(aGenoTable)
grdf <- as.data.frame(gr)
grdf

#Extracting the taxa(sample) wrapper from a GenotypeTable
aTaxa <- taxa(aGenoTable)
aTaxa
#Extract dataframe of taxa from GenotypeTable, then convert to data.frame
sampleDF <- sampleDataFrame(aGenoTable)
sampleDF

#create a summarizedExperiment  from a TASSEL read of the VCF
#SummarizedExperiment is a core data structure of Bioconductor
sumExp <- summarizeExperimentFromGenotypeTable(readGenotypeTable(vcfPath))
show(sumExp)
#create a snpMatrix of the snpStats from a TASSEL read of the VCF
asnpmat <- snpMatrixFromGenotypeTable(readGenotypeTable(vcfPath))

#PCA Analysis
xxmat <- xxt(asnpmat, correct.for.missing=FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]
evals <- evv$values[1:5]
evals
plot(pcs[,1],pcs[,2])

testTaxaFilterGT <- filterTaxaBuilderPlugin(aGenoTable,0.3, includeTaxa=TRUE)
taxa(testTaxaFilterGT)

testTaxaFilterGT <- filterTaxaBuilderPlugin(aGenoTable,0.3, includeTaxa=TRUE, taxaList ="M0297:C05F2ACXX:5:250021042,A659:C08L7ACXX:6:250048004")

##Phenotype - Genotype GLM tests
## VCF file path exmaple...
phenotypePath <- paste0(
  getwd(),
  "/data/mdp_traits.txt"
)
genotypePath <- paste0(
  getwd(),
  "/data/mdp_genotype.hmp.txt.gz"
)
gwasGeno <- readGenotypeTable(genotypePath)

#Need to convert to Pull function
jtsPhenotype <- new(J("net/maizegenetics/phenotype/PhenotypeBuilder"))$fromFile(phenotypePath)$build()$get(0L)

#Estimates BLUEs - not working anymore at boolean passing messed UP
blueReports <- fixedEffectLMPlugin(jtsPhenotype, phenoOnly=TRUE)

#Does GWAS after combining phenotype and genotype
genoPhenoCombined <- combineTasselGenotypePhenotype(gwasGeno@jtsGenotypeTable,jtsPhenotype)
gwasReports <- fixedEffectLMPlugin(genoPhenoCombined)

#GWAS reports contains two dataframes - one with marker tests, other with allele effects.

