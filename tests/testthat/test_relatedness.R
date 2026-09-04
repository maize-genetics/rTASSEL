# === Tests for relatedness methods =================================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasGeno <- rtObjs$gt_hmp


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

























