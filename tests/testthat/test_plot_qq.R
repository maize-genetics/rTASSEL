# === Tests for QQ plotting =========================================

test_that("plotQQ works correctly", {
    ## Load data ----
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    phenoPathFast <- system.file(
        "extdata",
        "mdp_traits_nomissing.txt",
        package = "rTASSEL"
    )
    tasGenoPhenoFast <- rTASSEL::readGenotypePhenotype(
        genoPathOrObj = genoPathHMP,
        phenoPathDFOrObj = phenoPathFast
    )
    fastRep <- rTASSEL::assocModelFitter(
        tasGenoPhenoFast,
        . ~ .,
        fastAssociation = TRUE,
        fitMarkers = TRUE,
        maxP = 1
    )
    tasBLUE <- rTASSEL::assocModelFitter(
        tasGenoPhenoFast,
        . ~ .,
        fitMarkers = FALSE
    )

    ## Test for errors ----
    expect_error(plotQQ(mtcars))
    expect_error(plotQQ(tasBLUE))

    ## Test for proper data return (static) ----
    testPlt <- plotQQ(fastRep)
    expect_true(object = is(testPlt, "gg"))
    testPlt <- plotQQ(fastRep, overlay = FALSE)
    testPltBld <- ggplot2::ggplot_build(testPlt)
    expect_equal(length(testPltBld$layout$facet_params$rows), 0)
    expect_equal(length(testPltBld$layout$facet_params$cols), 0)

    ## Test for proper data return (interactive) ----
    testPlt <- plotQQ(fastRep, interactive = TRUE)
    expect_true(is(testPlt, "plotly"))
    expect_true(is(testPlt, "htmlwidget"))
})


