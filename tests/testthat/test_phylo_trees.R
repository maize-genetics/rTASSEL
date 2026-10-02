# === Tests for phylogenetic analyses ===============================

skip_if_not_installed("ape")

## Preamble - load data ----

### Start logging info
startLogger()


### Shared fixtures (see helper_vars.R)
tasGeno    <- rtObjs$gt_hmp
tasGenoSub <- tasGeno[taxa("33-16", "38-11", "4226", "4722", "A188", "A214N"), ]


## Error tests ----
test_that("createTree() throws general exceptions.", {
    expect_error(
        object = createTree(tasObj = mtcars),
        regexp = "Unsupported input object"
    )
    expect_error(
        object = createTree(tasObj = tasGeno, clustMethod = "NA"),
        regexp = NULL
    )
})


## Return tests ----
test_that("createTree() returns correct data types", {
    t <- createTree(tasObj = tasGeno)
    tSub <- createTree(tasObj = tasGenoSub)

    expect_s3_class(
        object = t,
        class = "phylo"
    )

    expect_equal(
        object = t$Nnode,
        expected = 279
    )

    expect_equal(
        object = tSub$Nnode,
        expected = 4
    )

    expect_equal(
        object = length(t$tip.label),
        expected = 281
    )

    expect_equal(
        object = length(tSub$tip.label),
        expected = 6
    )
})


## Back-compatibility ----
test_that("createTree() accepts a deprecated TasselGenotypePhenotype", {
    expect_equal(
        createTree(tasObj = rtObjsLegacy$gt_hmp),
        createTree(tasObj = tasGeno)
    )
})


