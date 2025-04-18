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
rtObjs <- list(
    "gt_hmp"           = readGenotypeTableFromPath(rtFiles$gt_hmp_path),
    "gt_vcf"           = readGenotypeTableFromPath(rtFiles$gt_vcf_path),
    "ph_full"          = readPhenotypeFromPath(rtFiles$ph_full_path),
    "ph_nomiss"        = readPhenotypeFromPath(rtFiles$ph_nomiss_path),
    "gt_hmp_ph_full"   = readGenotypePhenotype(rtFiles$gt_hmp_path, rtFiles$ph_full_path),
    "gt_hmp_ph_nomiss" = readGenotypePhenotype(rtFiles$gt_hmp_path, rtFiles$ph_nomiss_path)
)


