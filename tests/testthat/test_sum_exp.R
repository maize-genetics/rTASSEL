# === Tests SummarizedExperiment creation method ====================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasPheno   <- rtObjs$ph_nomiss
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss


## Error tests ----
test_that("getSumExpFromGenotypeTable() throws general exceptions.", {
    expect_error(
        object = getSumExpFromGenotypeTable(
            tasObj            = mtcars,
            coerceDosageToInt = FALSE,
            verbose           = FALSE
        ),
        regexp = "Unsupported input object"
    )
    expect_error(
        object = getSumExpFromGenotypeTable(
            tasObj            = tasPheno,
            coerceDosageToInt = FALSE,
            verbose           = FALSE
        ),
        regexp = "needs genotype data"
    )
})


## Return tests ----
test_that("getSumExpFromGenotypeTable() returns correct data.", {
    tasSE <- getSumExpFromGenotypeTable(
        tasObj            = tasDataset,
        coerceDosageToInt = FALSE,
        verbose           = FALSE
    )

    expect_equal(
        object = class(tasSE)[1],
        expected = "RangedSummarizedExperiment"
    )
    expect_equal(
        object = dim(tasSE),
        expected = c(3093, 278)
    )
})

test_that("getSumExpFromGenotypeTable() returns correct dosage types.", {
    tasSERaw <- getSumExpFromGenotypeTable(
        tasObj            = tasDataset,
        coerceDosageToInt = FALSE,
        verbose           = FALSE
    )
    assaySERaw <- SummarizedExperiment::assay(tasSERaw)[[1, 1]]

    tasSEInt <- getSumExpFromGenotypeTable(
        tasObj            = tasDataset,
        coerceDosageToInt = TRUE,
        verbose           = FALSE
    )
    assaySEInt <- SummarizedExperiment::assay(tasSEInt)[[1, 1]]

    expect_equal(
        object = class(assaySEInt),
        expected = "integer"
    )
    expect_equal(
        object = class(assaySERaw),
        expected = "raw"
    )
})


## Back-compatibility ----
test_that("getSumExpFromGenotypeTable() accepts a deprecated TasselGenotypePhenotype", {
    tasSE <- getSumExpFromGenotypeTable(
        tasObj            = rtObjsLegacy$gt_hmp_ph_nomiss,
        coerceDosageToInt = FALSE,
        verbose           = FALSE
    )

    expect_equal(dim(tasSE), c(3093, 278))
})


