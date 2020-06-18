# === Tests for relatedness functions ===============================

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
tasGenoPheno <- readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoPathFast
)



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
        regexp = "`tasObj` must be of class `TasselGenotypePhenotype`"
    )
    expect_error(
        object = kinshipMatrix(
            tasObj             = tasPheno,
            method             = "Centered_IBS",
            maxAlleles         = 6,
            algorithmVariation = "Observed_Allele_Freq"
        ),
        regexp = "TASSEL genotype object not found"
    )
})

### distanceMatrix ----
test_that("distanceMatrix() throws general exceptions.", {
    expect_error(
        object = distanceMatrix(
            tasObj = mtcars
        ),
        regexp = "`tasObj` must be of class `TasselGenotypePhenotype`"
    )
    expect_error(
        object = distanceMatrix(
            tasObj = tasPheno
        ),
        regexp = "TASSEL genotype object not found"
    )
})



## Return tests ----

### kinshipMatrix ----
test_that("kinshipMatrix() returns correct data.", {
    tasKin <- kinshipMatrix(
        tasObj             = tasGenoPheno,
        method             = "Centered_IBS",
        maxAlleles         = 6,
        algorithmVariation = "Observed_Allele_Freq"
    )

    expect_equal(
        object = class(tasKin)[1],
        expected = "jobjRef"
    )
    expect_equal(
        object = length(names(tasKin)),
        expected = 35
    )


    tasKinRMat <- kinshipToRMatrix(tasKin)

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
    tasDist <- distanceMatrix(tasGenoPheno)

    expect_equal(
        object = class(tasDist)[1],
        expected = "jobjRef"
    )
    expect_equal(
        object = length(names(tasDist)),
        expected = 35
    )


    tasDistRMat <- distanceToRMatrix(tasDist)

    expect_equal(
        object   = class(tasDistRMat),
        expected = c("matrix", "array")
    )
    expect_equal(
        object   = dim(tasDistRMat),
        expected = c(278, 278)
    )

})


