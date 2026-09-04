# === Tests for `genomicPrediction()` ===============================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss

### Create kinship object
tasKin <- kinshipMatrix(tasDataset)



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
        regexp = "Unsupported input object"
    )
    expect_error(
        object = genomicPrediction(
            tasPhenoObj = tasGeno,
            kinship     = tasKin,
            doCV        = TRUE,
            kFolds      = 10,
            nIter       = 10
        ),
        regexp = "needs phenotype data"
    )
    expect_error(
        object = genomicPrediction(
            tasPhenoObj = tasDataset,
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
        tasPhenoObj = tasDataset,
        kinship     = tasKin,
        doCV        = TRUE,
        kFolds      = 2,
        nIter       = 1
    )
    gp <- genomicPrediction(
        tasPhenoObj = tasDataset,
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


## Back-compatibility ----
test_that("genomicPrediction() accepts a deprecated TasselGenotypePhenotype", {
    expect_equal(
        genomicPrediction(
            tasPhenoObj = rtObjsLegacy$gt_hmp_ph_nomiss,
            kinship     = tasKin,
            doCV        = FALSE
        ),
        genomicPrediction(
            tasPhenoObj = tasDataset,
            kinship     = tasKin,
            doCV        = FALSE
        )
    )
})











