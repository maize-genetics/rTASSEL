# Start rJava and set paths
setwd("/Users/sej65/rtassel")
path_tassel <- paste0(getwd(),"/inst/java/sTASSEL.jar")
Sys.setenv(JAVA_HOME="/Users/sej65/anaconda/")
library(rJava)

## jinit
rJava::.jinit(parameters="-Xmx12g")
.jcall(.jnew("java/lang/Runtime"), "J", "totalMemory")
.jcall(.jnew("java/lang/Runtime"), "J", "maxMemory")

## Set class paths
path_tassel <- "/Users/sej65/rtassel/inst/java/sTASSEL.jar"
path_phg <- "/Users/sej65/rtassel/inst/java/phg.jar"
path_kotlin <- "/Users/sej65/rtassel/inst/java/kotlin-stdlib-1.3.10.jar"
path_lib <- "/Users/sej65/rtassel/inst/java/lib/"
rJava::.jinit()
rJava::.jaddClassPath(path_tassel)
rJava::.jaddClassPath(path_phg)
rJava::.jaddClassPath(path_kotlin)
rJava::.jaddClassPath(path_lib)
print(.jclassPath())

# Same function as in PHGPluginWrappers.R
haplotypeGraphBuilderPlugin <- function(configFile, myMethods, chrom_list=NULL
) {
  plugin <- new(J("net.maizegenetics.pangenome.api.HaplotypeGraphBuilderPlugin"), .jnull(), FALSE)
  plugin$setParameter("configFile",toString(configFile))
  plugin$setParameter("methods",toString(myMethods))
  plugin$setParameter("includeSequences",toString(FALSE))
  plugin$setParameter("includeVariantContexts",toString(TRUE))
  plugin$setParameter("chromosomes", toString(chrom_list)) #comma separated list of chromosomes, taken as a string
  plugin$build()
}

#Create graph
configFilePath <- paste0(getwd(),"/data/configSQLiteR.txt")
method <- "mummer4,refRegionGroup"
test_graph <- haplotypeGraphBuilderPlugin(configFilePath,method)

# Call variant info
rrIds2 <- c(1:2)
generateRforPHG <- .jnew("net.maizegenetics.pangenome/pipelineTests.GenerateRForPHG")
hapDataVecs <- rJava::.jcall(generateRforPHG, "Lnet/maizegenetics/pangenome/pipelineTests/HaplotypesDataVectors;", "graphToHapsInRefRangeVectors", test_graph, rrIds2, FALSE, TRUE)
hapDataVecs$variantInfo

#Get genotype table
taxalist <- c("B104_Assembly")
chrom_list <- c(1)
genotypeTable <- J("net.maizegenetics/pangenome/api/GraphToGenotypeTable")$genotypeTable(test_graph, 1L, NULL)
genotypeTable2 <- J("net.maizegenetics/pangenome/api/GraphToGenotypeTable")$genotypeTable(test_graph, 2L, NULL)

# Convert genotype table to snpMatrix via RangedSummarizedExperiment class
# SNP_mx = snpMatrixFromGenotypeTable(genotypeTable) --> doesn't work see error
SumExp = summarizeExperimentFromGenotypeTable(genotypeTable = genotypeTable)
SNP_mx = assay(SumExp, 1)

SumExp2 = summarizeExperimentFromGenotypeTable(genotypeTable = genotypeTable2)
SNP_mx2 = assay(SumExp2, 1)
