# === Tests for scree plotting ======================================

test_that("plotScree works correctly", {
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

    ## General tests ----
    testPlt    <- plotScree(pcaRes)
    testPltBld <- ggplot2::ggplot_build(testPlt)
    expect_true(is(testPlt, "gg"))
    expect_error(plotScree(pcaRes, nComp = 500))
    expect_error(plotScree(mtcars))
    expect_error(
        plotScree(pca(tasGeno, reportEigenvalues = F, reportEigenvectors = F))
    )

    ## Test geometry ----
    testPlt    <- plotScree(pcaRes)
    testPltBld <- ggplot2::ggplot_build(testPlt)
    expect_equal(nrow(testPltBld$data[[1]]), 10)
    testPlt    <- plotScree(pcaRes, nComp = 20)
    testPltBld <- ggplot2::ggplot_build(testPlt)
    expect_equal(nrow(testPltBld$data[[1]]), 20)

    ## Test interactivity ----
    testPlt <- plotScree(pcaRes, nComp = 5, interactive = TRUE)
    expect_true(is(testPlt, "plotly"))
})

