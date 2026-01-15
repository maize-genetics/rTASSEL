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

test_that("genomicRanges validates inputs correctly", {
    expect_error(genomicRanges(mtcars), "data.frame does not contain a TASSEL PositionList object")
    expect_error(genomicRanges(tasPheno), "TasselGenotypePhenotype does not contain a TASSEL PositionList object")

    # Test with valid input
    gr <- genomicRanges(tasGeno)
    expect_s4_class(gr, "GRanges")
    expect_true("tasselIndex" %in% names(GenomicRanges::mcols(gr)))
})

test_that("tableReportToDF handles different input types", {
    # Test with NULL input
    expect_error(tableReportToDF(NULL))
})

test_that("checkForValidColumns validates column names", {
    test_df <- data.frame(a = 1:3, b = 4:6)
    expect_error(
        checkForValidColumns(test_df, c("a", "c")),
        "'c' column not found in stats dataframe"
    )

    # Should not error with valid columns
    expect_silent(checkForValidColumns(test_df, c("a", "b")))
})


# === Tests for Java build result helper functions ==================

test_that("safeGetFirst() returns NULL for empty or null Java lists", {
    # Test with NULL input
    expect_null(safeGetFirst(NULL))

    # Test with Java null
    expect_null(safeGetFirst(rJava::.jnull()))

    # Test with empty Java ArrayList
    emptyList <- rJava::.jnew("java.util.ArrayList")
    expect_null(safeGetFirst(emptyList))
})

test_that("safeGetFirst() returns first element for non-empty Java lists", {
    # Test with non-empty Java ArrayList
    javaList <- rJava::.jnew("java.util.ArrayList")
    javaList$add("first_element")
    javaList$add("second_element")

    result <- safeGetFirst(javaList)
    expect_equal(result, "first_element")
})

test_that("formatTaxaForMessage() formats taxa correctly", {
    # Test with empty taxa
    expect_equal(
        formatTaxaForMessage(character(0), "Test taxa"),
        "Test taxa: (none found)"
    )

    # Test with single taxon
    expect_equal(
        formatTaxaForMessage("taxon1", "Test taxa"),
        "Test taxa: 'taxon1'"
    )

    # Test with multiple taxa
    taxa <- c("taxon1", "taxon2", "taxon3")
    expect_equal(
        formatTaxaForMessage(taxa, "Test taxa"),
        "Test taxa: 'taxon1', 'taxon2', 'taxon3'"
    )

    # Test with more taxa than maxShow
    manyTaxa <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7")
    result <- formatTaxaForMessage(manyTaxa, "Test taxa", maxShow = 3)
    expect_true(grepl("\\.\\.\\.$", result))
    expect_true(grepl("'t1', 't2', 't3'", result))
})

test_that("phenotypeHasNoTaxa() correctly identifies empty phenotypes", {
    # Test with NULL input
    expect_true(phenotypeHasNoTaxa(NULL))

    # Test with Java null
    expect_true(phenotypeHasNoTaxa(rJava::.jnull()))

    # Test with valid phenotype (should have taxa)
    phenoPath <- system.file("extdata", "mdp_traits_nomissing.txt", package = "rTASSEL")
    validPheno <- readPhenotypeFromPath(phenoPath)
    expect_false(phenotypeHasNoTaxa(validPheno@jPhenotypeTable))
})

test_that("collectTaxaSamplesFromObjects() extracts taxa from objects", {
    # Create test phenotype objects
    phA <- readPhenotypeFromDataFrame(
        data.frame(taxa = c("a", "b", "c"), value = c(1, 2, 3)),
        "taxa"
    )
    phB <- readPhenotypeFromDataFrame(
        data.frame(taxa = c("x", "y"), value = c(10, 20)),
        "taxa"
    )

    taxaSamples <- collectTaxaSamplesFromObjects(list(phA, phB))

    expect_equal(length(taxaSamples), 2)
    expect_equal(taxaSamples[[1]], c("a", "b", "c"))
    expect_equal(taxaSamples[[2]], c("x", "y"))
})

test_that("extractTaxaSample() handles edge cases",
{
    # Test with NULL input
    expect_equal(extractTaxaSample(NULL), character(0))

    # Test with Java null
    expect_equal(extractTaxaSample(rJava::.jnull()), character(0))
})

test_that("extractTaxaSample() extracts taxa from valid TaxaList", {
    # Use existing genotype data to get a real TaxaList
    taxa <- extractTaxaSample(tasGeno@jTaxaList, maxSamples = 3)

    expect_type(taxa, "character")
    expect_equal(length(taxa), 3)
})


