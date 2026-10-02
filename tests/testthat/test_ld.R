# === Tests for `linkageDiseq()` ====================================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasPheno   <- rtObjs$ph_nomiss
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss



## Error tests ----

### General errors ----
test_that("linkageDiseq() throws general exceptions.", {
    expect_error(
        object = linkageDiseq(
            tasObj     = mtcars,
            ldType     = "slidingWindow",
            windowSize = NULL,
            hetCalls   = "missing"
        ),
        regexp = "Unsupported input object"
    )
    expect_error(
        object = linkageDiseq(
            tasObj     = tasPheno,
            ldType     = "slidingWindow",
            windowSize = NULL,
            hetCalls   = "missing"
        ),
        regexp = "needs genotype data"
    )
    expect_that(
        object = linkageDiseq(
            tasObj     = tasDataset,
            ldType     = "All",
            windowSize = NULL,
            hetCalls   = "mising"
        ),
        condition = throws_error()
    )
    expect_that(
        object = linkageDiseq(
            tasObj     = tasDataset,
            ldType     = "Everything",
            windowSize = NULL,
            hetCalls   = "missing"
        ),
        condition = throws_error()
    )

})



## Return tests ----
test_that("linkageDiseq() returns LDResults with correct structure.", {
    ldRes <- linkageDiseq(
        tasObj     = tasGeno,
        ldType     = "SlidingWindow",
        windowSize = 10,
        hetCalls   = "missing"
    )

    expect_s4_class(ldRes, "LDResults")
    expect_equal(ldRes@ldType, "SlidingWindow")
    expect_equal(ldRes@windowSize, 10)
    expect_equal(ldRes@hetCalls, "missing")

    ldDF <- ldRes@results
    expect_equal(dim(ldDF), c(30875, 17))

    expect_equal(
        object   = colnames(ldDF),
        expected = c(
            "Locus1", "Position1", "Site1", "NumberOfStates1", "States1",
            "Frequency1", "Locus2", "Position2", "Site2", "NumberOfStates2",
            "States2", "Frequency2", "Dist_bp", "R^2", "DPrime", "pDiseq", "N"
        )
    )
})

test_that("linkageDiseq() tableReport returns tibble.", {
    ldRes <- linkageDiseq(
        tasObj     = tasGeno,
        ldType     = "SlidingWindow",
        windowSize = 10,
        hetCalls   = "missing"
    )

    tbl <- tableReport(ldRes)
    expect_s3_class(tbl, "tbl_df")
    expect_equal(nrow(tbl), 30875)
})


## Back-compatibility ----
test_that("linkageDiseq() accepts a deprecated TasselGenotypePhenotype", {
    ldArgs <- list(
        ldType     = "SlidingWindow",
        windowSize = 10,
        hetCalls   = "missing"
    )

    legacyRes <- do.call(
        linkageDiseq,
        c(list(tasObj = rtObjsLegacy$gt_hmp), ldArgs)
    )
    modernRes <- do.call(
        linkageDiseq,
        c(list(tasObj = tasGeno), ldArgs)
    )

    expect_s4_class(legacyRes, "LDResults")
    expect_equal(legacyRes@results, modernRes@results)
})


