

test_that("pca() returns correct data", {
    ## Load data ----
    tasGeno <- rtObjs$gt_hmp

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


test_that("pca() accepts a genomic dataset", {
    pcaRes <- pca(rtObjs$ds_hmp_ph_nomiss)

    expect_true(is(pcaRes, "PCAResults"))
    expect_equal(
        nrow(tableReport(pcaRes, "PC_Datum")),
        length(taxaList(rtObjs$ds_hmp_ph_nomiss))
    )
})


## Back-compatibility ----
test_that("pca() accepts a deprecated TasselGenotypePhenotype", {
    expect_equal(
        tableReport(pca(rtObjsLegacy$gt_hmp), "PC_Datum"),
        tableReport(pca(rtObjs$gt_hmp), "PC_Datum")
    )
})






