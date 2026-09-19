# === Tests for MDSResults objects ==================================

test_that("MDSResults methods work correctly.", {
    ## Make test data ----
    set.seed(123)
    mockMDSResults <- list(
        "MDS_PCs_Datum" = data.frame(
            "Taxa" = letters[1:3],
            "PC1"  = rnorm(3),
            "PC2"  = rnorm(3),
            "PC3"  = rnorm(3)
        ),
        "MDS_Eigenvalues_Datum" = data.frame(
            "PC"         = c(1L, 2L, 3L),
            "eigenvalue" = c(1.37, 0.73, 0.37)
        )
    )

    ## Test instantiation ----
    testClass <- methods::new(
        "MDSResults",
        results = mockMDSResults
    )
    expect_true(is(testClass, "MDSResults"))

    ## Test validity ----
    expect_error(
        methods::new(
            "MDSResults",
            results = c(mockMDSResults, list("bad" = 1:3))
        ),
        regexp = "Invalid list"
    )

    ## Test show method ----
    testCapture <- capture.output(testClass)
    expect_true(any(grepl("MDSResults object with 2", testCapture)))
    expect_true(any(grepl("3 reported axes", testCapture)))
    expect_true(any(grepl("Results:", testCapture)))
    expect_true(any(grepl("MDS_PCs_Datum", testCapture)))
    expect_true(any(grepl("MDS_Eigenvalues_Datum", testCapture)))

    ## Test getter methods ----
    expect_equal(
        object = reportNames(testClass),
        expected = c("MDS_PCs_Datum", "MDS_Eigenvalues_Datum")
    )
    expect_equal(
        object = colnames(tableReport(testClass)),
        expected = c("Taxa", "PC1", "PC2", "PC3")
    )
    expect_equal(
        object = colnames(tableReport(testClass, "MDS_PCs_Datum")),
        expected = c("Taxa", "PC1", "PC2", "PC3")
    )
    expect_equal(
        object = colnames(tableReport(testClass, "MDS_Eigenvalues_Datum")),
        expected = c("PC", "eigenvalue")
    )
    expect_true(is.list(tableReport(testClass, "all")))
    expect_true(is.list(tableReport(testClass, "ALL")))
    expect_error(tableReport(testClass, "pca"))
    expect_error(tableReport(testClass, 123))
})
