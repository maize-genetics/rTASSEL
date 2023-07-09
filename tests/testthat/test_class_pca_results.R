# === Tests for PCAResults objects ==================================

test_that("PCAResults methods work correctly.", {
    ## Make test data ----
    set.seed(123)
    mockPCAResults <- list(
        "PC_Datum" = data.frame(
            "Taxa" = letters[1:3],
            "PC1"  = rnorm(3),
            "PC2"  = rnorm(3),
            "PC3"  = rnorm(3)
        ),
        "Eigenvalues_Datum" = data.frame(
            "PC"                    = c("0", "1", "2"),
            "eigenvalue"            = rnorm(3),
            "proportion_of_total"   = c(0.6, 0.3, 0.1),
            "cumulative_proportion" = c(0.6, 0.9, 1.0)
        ),
        "Eigenvectors_Datum" = data.frame(
            "trait"        = paste0("snp_", letters),
            "Eigenvector1" = rnorm(length(letters)),
            "Eigenvector2" = rnorm(length(letters)),
            "Eigenvector3" = rnorm(length(letters))
        )
    )

    ## Test instantiation ----
    testClass <- methods::new(
        "PCAResults",
        results = mockPCAResults
    )
    expect_true(is(testClass, "PCAResults"))

    ## Test show method ----
    testCapture <- capture.output(testClass)
    expect_true(any(grepl("PCAResults object with 3", testCapture)))
    expect_true(any(grepl("Results:", testCapture)))
    expect_true(any(grepl("PC_Datum", testCapture)))
    expect_true(any(grepl("Eigenvalues_Datum", testCapture)))

    ## Test getter methods ----
    expect_equal(
        object = reportNames(testClass),
        expected =  c("PC_Datum", "Eigenvalues_Datum", "Eigenvectors_Datum")
    )
    expect_equal(
        object = colnames(tableReport(testClass)),
        expected = c("Taxa", "PC1", "PC2", "PC3")
    )
    expect_equal(
        object = colnames(tableReport(testClass, "PC_Datum")),
        expected = c("Taxa", "PC1", "PC2", "PC3")
    )
    expect_true(is.list(tableReport(testClass, "all")))
    expect_true(is.list(tableReport(testClass, "ALL")))
    expect_error(tableReport(testClass, "mds"))
    expect_error(tableReport(testClass, 123))
})


