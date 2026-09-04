# === Tests for relatedness functions ===============================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasPheno   <- rtObjs$ph_nomiss
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss



## Error tests ----

### kinshipMatrix ----
test_that("kinshipMatrix() throws general exceptions.", {
    expect_error(
        object = kinshipMatrix(
            tasObj             = mtcars,
            method             = "Centered_IBS",
            maxAlleles         = 6,
            algorithmVariation = "Observed_Allele_Freq"
        ),
        regexp = "Unsupported input object"
    )
    expect_error(
        object = kinshipMatrix(
            tasObj             = tasPheno,
            method             = "Centered_IBS",
            maxAlleles         = 6,
            algorithmVariation = "Observed_Allele_Freq"
        ),
        regexp = "needs genotype data"
    )
})

### distanceMatrix ----
test_that("distanceMatrix() throws general exceptions.", {
    expect_error(
        object = distanceMatrix(
            tasObj = mtcars
        ),
        regexp = "Unsupported input object"
    )
    expect_error(
        object = distanceMatrix(
            tasObj = tasPheno
        ),
        regexp = "needs genotype data"
    )
})



## Return tests ----

### kinshipMatrix ----
test_that("kinshipMatrix() returns correct data.", {
    tasKin <- kinshipMatrix(
        tasObj             = tasDataset,
        method             = "Centered_IBS",
        maxAlleles         = 6,
        algorithmVariation = "Observed_Allele_Freq"
    )

    expect_equal(
        object = class(tasKin)[1],
        expected = "TasselDistanceMatrix"
    )


    tasKinRMat <- as.matrix(tasKin)

    expect_equal(
        object   = class(tasKinRMat),
        expected = c("matrix", "array")
    )
    expect_equal(
        object   = dim(tasKinRMat),
        expected = c(278, 278)
    )

})


### distanceMatrix ----
test_that("distanceMatrix() returns correct data.", {
    tasDist <- distanceMatrix(tasDataset)

    expect_equal(
        object = class(tasDist)[1],
        expected = "TasselDistanceMatrix"
    )

    tasDistRMat <- as.matrix(tasDist)

    expect_equal(
        object   = class(tasDistRMat),
        expected = c("matrix", "array")
    )
    expect_equal(
        object   = dim(tasDistRMat),
        expected = c(278, 278)
    )

})


### Bare genotype input ----
test_that("relatedness methods accept a genotype without phenotype data.", {
    expect_equal(dim(as.matrix(kinshipMatrix(tasGeno))), c(281, 281))
    expect_equal(dim(as.matrix(distanceMatrix(tasGeno))), c(281, 281))
})


## Back-compatibility ----
test_that("relatedness methods accept a deprecated TasselGenotypePhenotype", {
    legacyObj <- rtObjsLegacy$gt_hmp_ph_nomiss

    expect_equal(
        as.matrix(kinshipMatrix(legacyObj)),
        as.matrix(kinshipMatrix(tasDataset))
    )
    expect_equal(
        as.matrix(distanceMatrix(legacyObj)),
        as.matrix(distanceMatrix(tasDataset))
    )
})


