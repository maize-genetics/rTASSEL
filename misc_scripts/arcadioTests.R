# Test script for creating SummarizedExperimemt as per
# https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html

# extends original initialTests.R and btmTests.R of converting tassel GenotypeTable to GenomicRanges

library(rJava)
library(GenomicRanges)
library(SummarizedExperiment)
library(GWASpoly)
library(ggplot2)

is_experimental <- TRUE

#set workdir for rtassel
setwd("~/myBins/bucklerlabBitbucket/rtassel/")

path_exp_tassel <- paste0(getwd(),"/inst/java/sTASSEL.jar")
path_exp_tassel_libs <- paste0(getwd(),"/inst/java/lib")

## jinit
rJava::.jinit(parameters="-Xmx6g")
.jcall(.jnew("java/lang/Runtime"), "J", "totalMemory")
.jcall(.jnew("java/lang/Runtime"), "J", "maxMemory")

## Add class paths
if(is_experimental == TRUE) {
  tasselPath <- path_exp_tassel
  tasselLibs <- path_exp_tassel_libs
}

rJava::.jaddClassPath(tasselPath)
rJava::.jaddClassPath(tasselLibs)
print(.jclassPath())


## Source files
source("R/AllGenerics.R")
source("R/AllClasses.R")
source("R/TasselPluginWrappers.R")
source("R/PullFunctions.R")
source("R/gwasPolyObjectCreator.R")

## Add VCF path
vcfPath <- paste0(
  getwd(),
  "/data/maize_chr9_10thin40000.recode.vcf"
)


# Tests

## Make genotype table from tasses sample data
tasGenoTable <- readGenotypeTable(paste0(getwd(), "/data/mdp_genotype.hmp.txt"))

## Make summarized experiment from genotypetable
tas_se <- summarizeExperimentFromGenotypeTable(tasGenoTable) # not working right now, but

## Load phenotype data
###straight load as dataframe, skpping first two rows on tassel specific phenotype table format
phenos <- read.table(file = paste0(getwd(), "/data/mdp_phenotype.txt"), skip = 2, header = T, sep = "\t", na.strings = "NaN")

### select sinlge location, as GWASpoly requires single entries for taxa.
phenosOneLoc <- phenos[phenos$location == "A",]
rownames(phenosOneLoc) <- phenosOneLoc$Taxa
###remove location as it is now redundant. Also, GWASpoly expects all traits as initial columns, and fixed effect covariates last
phenosOneLoc <- phenosOneLoc[,-c(2)]

## create GWASpoly object with coopted read.GWASpoly function
data_gwasPoly <- se_createGWASpolyObject(ploidy = 2, phenoDF = phenosOneLoc, SummarizedExperimentObject =  tas_se, format = "numeric", n.traits = 3)

## add kinship information to object
data_gwasPoly <- set.K(data_gwasPoly)

## set parameters for mixed model
params <- set.params(fixed=unlist(strsplit("Q1,Q2,Q3", ",")),
                     fixed.type=rep("numeric",3))

## run gwas
data_gwasPoly_res <- GWASpoly(data = data_gwasPoly, models = "additive", params = params)

qq.plot(data_gwasPoly_res, trait = "EarDia", model = "additive")

data_gwasPoly_res <- set.threshold(data_gwasPoly_res, method = "FDR", level=0.05)

get.QTL(data = data_gwasPoly_res)

manhattan.plot(data_gwasPoly_res, trait = "EarDia", model = "additive")

traitGWASresults <- data.frame()
for(trait in names(data_gwasPoly_res@scores)){
  message(trait)
  traitMarkerpScores <- as.data.frame(data_gwasPoly_res@scores[trait])
  Marker <- rownames(traitMarkerpScores)
  colnames(traitMarkerpScores) <- "markerLogPVal"
  traitMarkerpVals <- exp(traitMarkerpScores*-1)
  colnames(traitMarkerpVals) <- "markerpVal"
  traitMarkerEffects <-  as.data.frame(data_gwasPoly_res@effects[trait])
  colnames(traitMarkerEffects) <- "markerEffect"
  sigTreshold <- data_gwasPoly_res@threshold[rownames(data_gwasPoly_res@threshold)==trait]
  traitGWASresults <- rbind(traitGWASresults, data.frame(Marker, traitMarkerEffects, traitMarkerpVals, traitMarkerpScores, trait, sigTreshold))
}

summary(traitGWASresults)

traitGWASresults <- merge(traitGWASresults, data_gwasPoly_res@map, by = "Marker")
summary(traitGWASresults)

ggplot(traitGWASresults) + geom_point(aes(Position, markerLogPVal, col = trait), alpha=0.3) + geom_point(aes(Position, markerLogPVal, col = trait), data = traitGWASresults[traitGWASresults$markerLogPVal>traitGWASresults$sigTreshold,])  + geom_hline(aes(yintercept = sigTreshold, col = trait), lty = "dashed") + facet_grid(~Chrom, scales = "free_x") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_brewer(palette = "Set2") + xlab("Position by Chr") + ylab("-Log10Pval")

