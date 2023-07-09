# === Tests for relatedness methods =================================

## Preamble - load data ----

### Start logging info
startLogger()

### Load hapmap data
genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)

### Read data - need only non missing data!
phenoPathFast <- system.file(
    "extdata",
    "mdp_traits_nomissing.txt",
    package = "rTASSEL"
)

### Create rTASSEL phenotype only object
tasPheno <- readPhenotypeFromPath(
    path = phenoPathFast
)

### Create rTASSEL genotype only object
tasGeno <- readGenotypeTableFromPath(
    path = genoPathHMP
)

### Create rTASSEL object - use prior TASSEL genotype object
tasGenoPhenoFast <- readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoPathFast
)


test_that("asTasselDistanceMatrix() returns correct data and exceptions", {
    set.seed(123)
    m <- 10
    s <- matrix(rnorm(100), m)
    s[lower.tri(s)] <- t(s)[lower.tri(s)]
    diag(s) <- 2

    expect_error(
        object = asTasselDistanceMatrix(s),
        regexp = "Matrix object must have column and row"
    )
})


test_that("mds() returns correct data and exceptions", {
    expect_error(mds(tasGeno))

    # set.seed(123)
    # m <- 10
    # s <- matrix(rnorm(100), m)
    # s[lower.tri(s)] <- t(s)[lower.tri(s)]
    # diag(s) <- 2
    # colnames(s) <- rownames(s) <- paste0("s_", seq_len(m))
    #
    # sT <- distanceMatrix(tasGeno)
    #
    # expect_true(inherits(mds(sT, "list")))
})

























