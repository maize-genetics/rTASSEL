# === Tests for `genomicPrediction()` ===============================

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

### Create rTASSEL genotype only object
tasGeno <- readGenotypeTableFromPath(
    path = genoPathHMP
)

### Create rTASSEL object - use prior TASSEL genotype object
tasGenoPheno <- readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoPathFast
)

### Create kinship object
tasKin <- kinshipMatrix(tasGenoPheno)



## Error tests ----
test_that("genomicPrediction() throws general exceptions.", {
    expect_error(
        object = genomicPrediction(
            tasPhenoObj = mtcars,
            kinship     = tasKin,
            doCV        = TRUE,
            kFolds      = 10,
            nIter       = 10
        ),
        regexp = "`tasObj` must be of class `TasselGenotypePhenotype`"
    )
    expect_error(
        object = genomicPrediction(
            tasPhenoObj = tasGeno,
            kinship     = tasKin,
            doCV        = TRUE,
            kFolds      = 10,
            nIter       = 10
        ),
        regexp = "TASSEL phenotype object not found"
    )
    expect_error(
        object = genomicPrediction(
            tasPhenoObj = tasGenoPheno,
            kinship     = mtcars,
            doCV        = TRUE,
            kFolds      = 10,
            nIter       = 10
        ),
        regexp = "TASSEL kinship object is not of TasselDistanceMatrix class"
    )
})



## Return tests ----
test_that("genomicPrediction() returns correct data.", {
    gpCV <- genomicPrediction(
        tasPhenoObj = tasGenoPheno,
        kinship     = tasKin,
        doCV        = TRUE,
        kFolds      = 2,
        nIter       = 1
    )
    gp <- genomicPrediction(
        tasPhenoObj = tasGenoPheno,
        kinship     = tasKin,
        doCV        = FALSE
    )

    expect_equal(
        object   = colnames(gpCV),
        expected = c("Trait", "Iteration", "Fold", "Accuracy")
    )
    expect_equal(
        object   = colnames(gp),
        expected = c("Trait", "Taxon", "Observed", "Predicted", "PEV")
    )
    expect_equal(
        object = dim(gpCV),
        expected = c(6, 4)
    )
    expect_equal(
        object = dim(gp),
        expected = c(834, 5)
    )
})











