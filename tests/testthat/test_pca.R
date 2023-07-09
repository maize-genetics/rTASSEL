

test_that("pca() returns correct data", {
    ## Load data ----
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGeno <- readGenotypeTableFromPath(
        path = genoPathHMP
    )

    ## Run test PCA ----
    pcaRes <- pca(tasGeno)

    ## Tests ----
    expect_true(is(pcaRes, "PCAResults"))
    expect_equal(
        reportNames(pcaRes),
        c("PC_Datum", "Eigenvalues_Datum", "Eigenvectors_Datum")
    )

    expect_equal(
        reportNames(
            pca(
                tasGeno,
                reportEigenvalues = FALSE,
                reportEigenvectors = TRUE
            )
        ),
        c("PC_Datum", "Eigenvectors_Datum")
    )
    expect_equal(
        reportNames(
            pca(
                tasGeno,
                reportEigenvalues = TRUE,
                reportEigenvectors = FALSE
            )
        ),
        c("PC_Datum", "Eigenvalues_Datum")
    )
    expect_equal(
        reportNames(
            pca(
                tasGeno,
                reportEigenvalues = FALSE,
                reportEigenvectors = FALSE
            )
        ),
        "PC_Datum"
    )
    expect_equal(
        rJava::.jclass(pcaRes@jObj),
        "net.maizegenetics.phenotype.CorePhenotype"
    )
})






