# === Tests for utility methods =====================================
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
tasGenoPhenoFast <- readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoPathFast
)


test_that("truncate() function returns correct output", {
    s <- "This is a very long string!"
    sLenTruth <- nchar(s)

    expect_equal(truncate(s), "This is...")
    expect_equal(nchar(truncate(s)), 10)
    expect_equal(nchar(truncate(s, max = 25)), 25)
    expect_equal(truncate(s, max = 100), "This is a very long string!")
    expect_equal(nchar(truncate(s, max = 100)), sLenTruth)
    expect_equal(truncate(s, etc = "***"), "This is***")

    expect_error(truncate(s, max = -20))
})

test_that("summaryDistance() ", {
    ## Create a dummy pairwise matrix object ----
    set.seed(123)
    m <- 100
    s <- matrix(rnorm(m * m), m)
    s[lower.tri(s)] <- t(s)[lower.tri(s)]
    diag(s) <- 2

    ## Add sample IDs ----
    colnames(s) <- rownames(s) <- paste0("s_", seq_len(m))

    testTasselDist <- asTasselDistanceMatrix(s)

    expect_equal(ncol(summaryDistance(testTasselDist@jDistMatrix)), 7)
    expect_equal(nrow(summaryDistance(testTasselDist@jDistMatrix)), 7)

    ## Create a dummy pairwise matrix object ----
    set.seed(123)
    m <- 3
    s <- matrix(rnorm(m * m), m)
    s[lower.tri(s)] <- t(s)[lower.tri(s)]
    diag(s) <- 2

    ## Add sample IDs ----
    colnames(s) <- rownames(s) <- paste0("s_", seq_len(m))

    testTasselDist <- asTasselDistanceMatrix(s)

    expect_equal(ncol(summaryDistance(testTasselDist@jDistMatrix)), 4)
    expect_equal(nrow(summaryDistance(testTasselDist@jDistMatrix)), 4)

    ## Create a dummy pairwise matrix object ----
    set.seed(123)
    m <- 5
    s <- matrix(rnorm(m * m), m)
    s[lower.tri(s)] <- t(s)[lower.tri(s)]
    diag(s) <- 2

    ## Add sample IDs ----
    colnames(s) <- rownames(s) <- paste0("s_", seq_len(m))

    testTasselDist <- asTasselDistanceMatrix(s)

    expect_equal(ncol(summaryDistance(testTasselDist@jDistMatrix)), 6)
    expect_equal(nrow(summaryDistance(testTasselDist@jDistMatrix)), 6)
})

test_that("getGenotypePhenotype() returns correct data", {
    testObj <- getGenotypePhenotype(tasGenoPhenoFast)
    expect_equal(
        object = testObj$getClass()$toString(),
        expected = "class net.maizegenetics.phenotype.GenotypePhenotype"
    )

    testObj <- getGenotypePhenotype(tasGeno)
    expect_true(rJava::is.jnull(testObj))
})

test_that("assocStatsColumnChecker returns correct errors", {
    expect_error(
        object = tableReportListToAssociationResults(
            trl = list(),
            "FFF"
        ),
        regexp = "Association Type"
    )
})


