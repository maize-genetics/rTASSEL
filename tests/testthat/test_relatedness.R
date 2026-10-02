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

    tasDist <- distanceMatrix(tasGeno)
    mdsRes  <- mds(tasDist)

    expect_s4_class(mdsRes, "MDSResults")
    expect_equal(
        reportNames(mdsRes),
        c("MDS_PCs_Datum", "MDS_Eigenvalues_Datum")
    )
    expect_equal(
        colnames(tableReport(mdsRes)),
        c("Taxa", paste0("PC", 1:5))
    )
    expect_equal(
        nrow(tableReport(mdsRes)),
        length(taxaList(tasGeno))
    )
    expect_equal(
        nrow(tableReport(mdsRes, "MDS_Eigenvalues_Datum")),
        5
    )
    expect_equal(
        rJava::.jclass(mdsRes@jObj),
        "net.maizegenetics.phenotype.CorePhenotype"
    )

    expect_equal(
        colnames(tableReport(mds(tasDist, nAxes = 3))),
        c("Taxa", paste0("PC", 1:3))
    )
})

























