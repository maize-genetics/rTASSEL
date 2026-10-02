# === Tests for taxa lists ==========================================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss


## Tests ----
test_that("getTaxaList() returns correct data", {
    testObj <- getTaxaList(mtcars)
    expect_true(rJava::is.jnull(testObj))

    # A position list carries no taxa
    testObj <- getTaxaList(getPositionList(tasGeno))
    expect_true(rJava::is.jnull(testObj))
})

test_that("getTaxaList() accepts every class carrying taxa", {
    expect_equal(getTaxaList(tasGeno)$size(), 281L)
    expect_equal(getTaxaList(tasDataset)$size(), 278L)
    expect_equal(getTaxaList(rtObjs$ph_nomiss)$size(), 298L)
    expect_equal(getTaxaList(rtObjsLegacy$gt_hmp_ph_nomiss)$size(), 278L)
})

test_that("taxaList() returns correct excpetions", {
    expect_error(taxaList(mtcars))
})

test_that("getTaxaIDs() returns correct excpetions", {
    expect_error(getTaxaIDs(mtcars))
})

test_that("getTaxaIDs() agrees with taxaList() across classes", {
    expect_equal(getTaxaIDs(tasGeno), taxaList(tasGeno))
    expect_equal(getTaxaIDs(tasDataset), taxaList(tasDataset))
    expect_equal(getTaxaIDs(rtObjsLegacy$gt_hmp), taxaList(tasGeno))
})
