## ----
# TASSEL JVM Constants
# 
# @description
# Constants for TASSEL 5 Java class references used throughout the 
# package. These constants provide mappings to the Java classes and 
# methods used by rTASSEL.
TASSEL_JVM <- list(
    "ARRAY_LIST"               = "java.util.ArrayList",
    "CHROMOSOME"               = "net.maizegenetics.dna.map.Chromosome",
    "GENERAL_POSITION_BUILDER" = "net.maizegenetics.dna.map.GeneralPosition$Builder",
    "GENOTYPE_TABLE_BUILDER"   = "net.maizegenetics.dna.snp.GenotypeTableBuilder",
    "LOGGING_UTILS"            = "net.maizegenetics.util.LoggingUtils",
    "PHENO_BUILDER"            = "net.maizegenetics.phenotype.PhenotypeBuilder",
    "POSITION_LIST_BUILDER"    = "net.maizegenetics.dna.map.PositionListBuilder",
    "R_METHODS"                = "net.maizegenetics.plugindef.GenerateRCode",
    "REF_PROBABILITY_BUILDER"  = "net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder",
    "STEPWISE_FITTER"          = "net.maizegenetics.analysis.modelfitter.StepwiseOLSModelFitter",
    "STEPWISE_PLUGIN"          = "net.maizegenetics.analysis.modelfitter.StepwiseOLSModelFitterPlugin",
    "TAXA_LIST_BUILDER"        = "net.maizegenetics.taxa.TaxaListBuilder"
)


