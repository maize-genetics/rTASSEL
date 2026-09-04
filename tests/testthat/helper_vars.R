# /// Helper functions //////////////////////////////////////////////

## ----
# Return commonly used external package files
returnSysFiles <- function(x) system.file("extdata", x, package = "rTASSEL")


## ----
# Simulate numeric matrices
simNumericGt <- function(nRow, nCol) {
    # Matrix values
    minMax <- function(x) (x - min(x)) / (max(x) - min(x))
    d <- rnorm(nCol * nRow) |> minMax()
    m <- matrix(d, nrow = nRow, ncol = nCol)

    # Taxa values
    taxa <- sprintf("line_%02d", seq_len(nRow))

    # Position values
    mIds <- sprintf("marker_%02d", seq_len(nCol))
    mPos <- seq_len(nCol)

    # Add IDs to matrix
    colnames(m) <- mIds
    rownames(m) <- taxa

    return(m)
}



# /// Options ///////////////////////////////////////////////////////

## ----
# The 0.14.0 soft deprecations are asserted on purpose in
# test_deprecations.R, where `lifecycle::expect_deprecated()` re-enables
# them locally. Silencing them here keeps the rest of the suite - and the
# legacy fixtures below - free of notices that are not what is under test.
options(lifecycle_verbosity = "quiet")



# /// Constants /////////////////////////////////////////////////////

## ----
# Commonly used file paths
rtFiles <- list(
    "gt_hmp_path"       = returnSysFiles("mdp_genotype.hmp.txt"),
    "gt_vcf_path"       = returnSysFiles("maize_chr9_10thin40000.recode.vcf"),
    "ph_full_path"      = returnSysFiles("mdp_phenotype.txt"),
    "ph_nomiss_path"    = returnSysFiles("mdp_traits_nomissing.txt"),
    "ph_popstruct_path" = returnSysFiles("mdp_population_structure.txt")
)


## ----
# R matrices for numeric genotypes
rtMatrices <- list(
    "num_gt_sm" = simNumericGt(3, 3),
    "num_gt_md" = simNumericGt(10, 10),
    "num_gt_lg" = simNumericGt(50, 50)
)


## ----
# General rTASSEL objects
#
# The primary (0.14.0) classes: TasselGenotype, TasselPhenotype, and
# TasselGenomicDataset. The genotype and phenotype objects are reused when
# building the datasets so that the hapmap file is only parsed once.
rtObjs <- local({
    gtHmp    <- readGenotype(rtFiles$gt_hmp_path)
    phFull   <- readPhenotype(rtFiles$ph_full_path)
    phNoMiss <- readPhenotype(rtFiles$ph_nomiss_path)

    list(
        "gt_hmp"           = gtHmp,
        "gt_vcf"           = readGenotype(rtFiles$gt_vcf_path),
        "ph_full"          = phFull,
        "ph_nomiss"        = phNoMiss,
        "ds_hmp_ph_full"   = readGenomicDataset(gtHmp, phFull),
        "ds_hmp_ph_nomiss" = readGenomicDataset(gtHmp, phNoMiss)
    )
})


## ----
# Deprecated TasselGenotypePhenotype objects
#
# Kept so that each suite can hold on to a single back-compatibility case
# asserting that legacy pipelines still run and still get legacy classes back.
rtObjsLegacy <- local({
    gtHmp <- readGenotypeTableFromPath(rtFiles$gt_hmp_path)

    list(
        "gt_hmp"           = gtHmp,
        "ph_nomiss"        = readPhenotypeFromPath(rtFiles$ph_nomiss_path),
        "gt_hmp_ph_full"   = readGenotypePhenotype(gtHmp, rtFiles$ph_full_path),
        "gt_hmp_ph_nomiss" = readGenotypePhenotype(gtHmp, rtFiles$ph_nomiss_path)
    )
})


