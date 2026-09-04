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
    "GENO_PHENO_BUILDER"       = "net.maizegenetics.phenotype.GenotypePhenotypeBuilder",
    "GENOTYPE_TABLE_BUILDER"   = "net.maizegenetics.dna.snp.GenotypeTableBuilder",
    "LOGGING_UTILS"            = "net.maizegenetics.util.LoggingUtils",
    "PHENO_BUILDER"            = "net.maizegenetics.phenotype.PhenotypeBuilder",
    "POSITION_LIST_BUILDER"    = "net.maizegenetics.dna.map.PositionListBuilder",
    "R_METHODS"                = "net.maizegenetics.plugindef.GenerateRCode",
    "REF_PROBABILITY_BUILDER"  = "net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder",
    "STEPWISE_FITTER"          = "net.maizegenetics.analysis.modelfitter.StepwiseOLSModelFitter",
    "STEPWISE_PLUGIN"          = "net.maizegenetics.analysis.modelfitter.StepwiseOLSModelFitterPlugin",
    "TAXA_LIST_BUILDER"        = "net.maizegenetics.taxa.TaxaListBuilder",
    "VERSIONS"                 = "net.maizegenetics.tassel.TasselVersions"
)



## ----
# TASSEL Maven Jar Constants
#
# @description
# Constants for downloading the necessary JARs for internal TASSEL calls from
# Maven
TASSEL_MAVEN <- list(
    "BASE_URL"      = "https://repo1.maven.org/maven2",
    "GROUP_PATH"    = "net/maizegenetics",
    "GROUP_ID"      = "net.maizegenetics",
    "ARTIFACT_ID"   = "tassel",
    "VERSION"       = "5.2.96",
    "CLASSIFIER"    = "jar-with-dependencies",
    "SHA1_CHECKSUM" = "9320966721a12741da2a60f02fd3830639058d63",

    # TASSEL's own build declares SciJava (for 'cisd:jhdf5') and JBoss (for
    # 'openchart', a transitive of 'forester') in addition to Maven Central.
    # Neither artifact is resolvable from Central alone.
    "REPOS" = c(
        "central" = "https://repo1.maven.org/maven2",
        "scijava" = "https://maven.scijava.org/content/repositories/public",
        "jboss"   = "https://repository.jboss.org/maven2"
    )
)


## ----
# TASSEL GitHub Release Constants
#
# @description
# Constants for installing TASSEL from the standalone archives published on
# the TASSEL GitHub releases page. Nightly builds are only published there,
# so this is the only route to a bleeding-edge TASSEL.
#
# Each archive unpacks to a flat directory holding 'sTASSEL.jar' alongside a
# 'lib' directory of dependency JARs.
TASSEL_GITHUB <- list(
    "API_BASE"     = "https://api.github.com",
    "REPO"         = "maize-genetics/tassel",
    "ASSET_PREFIX" = "tassel-5-standalone-v",
    "ASSET_EXT"    = ".tar.gz",
    "MAIN_JAR"     = "sTASSEL.jar",
    "LIB_DIR"      = "lib",
    "PER_PAGE"     = 30L,
    "CACHE_FILE"   = "nightly-check.rds",
    "TIMEOUT_SECS" = 10L
)


## ----
# TASSEL Update Check Constants
#
# @description
# Settings governing the Maven Central version check performed when the
# package is attached.
TASSEL_UPDATE <- list(
    "CACHE_FILE"   = "version-check.rds",
    "MAX_AGE_SECS" = 24 * 60 * 60,
    "TIMEOUT_SECS" = 3L
)


