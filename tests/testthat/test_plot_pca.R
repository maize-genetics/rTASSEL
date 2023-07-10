# === Tests for PCA plotting ========================================

test_that("plotPCA works correctly.", {
    ## Load data ----
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    metadataPath <- system.file(
        "extdata",
        "mdp_metadata.csv",
        package = "rTASSEL"
    )
    tasGeno <- readGenotypeTableFromPath(
        path = genoPathHMP
    )
    tasMeta <- read.csv(metadataPath)
    tasMetaFaulty <- tasMeta
    colnames(tasMetaFaulty) <- c("sample", "Subpopulation")
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
    faultyPCARes <- methods::new(
        "PCAResults",
        results = mockPCAResults[c(F, T, T)]
    )

    ## Run test PCA ----
    pcaObj  <- pca(tasGeno)
    pcaObj2 <- pca(tasGeno, reportEigenvalues = FALSE)

    ## General tests ----
    testPlt    <- plotPCA(pcaObj)
    testPltBld <- ggplot2::ggplot_build(testPlt)
    expect_true(is(testPlt, "gg"))
    expect_true(is(plotPCA(pcaObj, interactive = TRUE), "plotly"))
    expect_error(plotPCA(mtcars))
    expect_error(plotPCA(faultyPCARes))
    expect_error(plotPCA(pcaObj, metadata = tasMeta, mCol = "error"))
    expect_error(plotPCA(pcaObj, 1, 20))
    expect_error(plotPCA(pcaObj, 1, "PC20"))
    expect_error(plotPCA(pcaObj, "PC100", 2))
    expect_error(plotPCA(pcaObj, "PCDF", 2))
    expect_error(plotPCA(pcaObj, metadata = iris, mCol = "Species"))
    expect_message(plotPCA(pcaObj2))
    expect_true(is(plotPCA(pcaObj, cluster = TRUE), "gg"))
    expect_true(is(plotPCA(pcaObj, metadata = tasMeta, mCol = "Subpopulation"), "gg"))
    expect_error(plotPCA(pcaObj, metadata = tasMetaFaulty, mCol = "Subpopulation"))
})





