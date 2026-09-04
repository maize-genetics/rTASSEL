# === Tests for position list functions =============================

## Preamble - load data ----

### Start logging info
startLogger()
library(rJava)

### Shared fixtures (see helper_vars.R)
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss


## Error tests ----
test_that("getPositionList() throws general exceptions.", {
    tmp <- getPositionList(mtcars)
    expect_true(rJava::is.jnull(tmp))
})


## Return tests ----
test_that("getPositionList() accepts every class carrying positions.", {
    isPositionList <- function(x) {
        getPositionList(x) %instanceof% "net.maizegenetics.dna.map.PositionArrayList"
    }

    expect_true(isPositionList(tasGeno))
    expect_true(isPositionList(tasDataset))
    expect_true(isPositionList(rtObjsLegacy$gt_hmp))

    # Raw Java references pass straight through
    expect_true(isPositionList(getPositionList(tasGeno)))
})

test_that("position list methods return correct data and classes.", {
    tmp <- genomicRanges(tasGeno)
    expect_true(class(tmp)[1] == "GRanges")
    expect_true(length(tmp$tasselIndex) == 3093)
})

test_that("positionList() agrees across the classes holding the table.", {
    expect_equal(positionList(tasGeno), positionList(tasDataset))
    expect_equal(nrow(positionList(tasGeno)), 3093)
})
