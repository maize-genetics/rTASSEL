# === Tests for Manhattan Plot functions ============================

## Equality tests ----
test_that("plotManhattan returns correct plot layers", {

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

    ## Check general errors ----
    expect_error(
        object = plotManhattan(mtcars),
        regexp = "The object 'mtcars' is not an 'AssociationResults' object"
    )
    expect_error(
        object = plotManhattan(tasBLUE),
        regexp = "Association Type not defined"
    )
    expect_error(
        object = plotManhattan(fastRep, "dog"),
        regexp = "No traits specified are found in results"
    )

    ## Check colors ----
    testPlt <- plotManhattan(
        assocRes = fastRep,
        trait = "dpoll",
        threshold = 2
    )
    testPltBld <- ggplot2::ggplot_build(testPlt)
    expect_true(is(testPlt, "gg"))
    expect_true(
        all(testPltBld$data[[1]]["colour"] |> unique() == c("#91baff", "#3e619b"))
    )
    testPlt <- plotManhattan(
        assocRes = fastRep,
        trait = "dpoll",
        threshold = 2,
        colors = c("#111111", "#444444")
    )
    testPltBld <- ggplot2::ggplot_build(testPlt)
    expect_true(is(testPlt, "gg"))
    expect_true(
        all(testPltBld$data[[1]]["colour"] |> unique() == c("#111111", "#444444"))
    )

    ## Test for ggplot object components ----
    testPlt <- plotManhattan(
        assocRes  = fastRep,
        trait     = "EarHT",
        threshold = 6
    )
    expect_s3_class(
        object = testPlt$layers[[1]],
        class  = "ggproto"
    )
    expect_s3_class(
        object = testPlt$layers[[1]]$geom,
        class  = "GeomPoint"
    )
    expect_equal(
        object   = testPlt$labels$title,
        expected = "Trait: EarHT"
    )
    expect_equal(
        object   = testPlt$labels$x,
        expected = "SNP Position (Mbp)"
    )
    expect_equal(
        object   = testPlt$labels$y,
        expected = ~-log[10] ~ "(" * italic(p) * "-value)"
    )

    ## Test for ggplot facets ----
    testPlt1 <- ggplot2::ggplot_build(plotManhattan(fastRep))
    testPlt2 <- ggplot2::ggplot_build(plotManhattan(fastRep, "dpoll"))
    expect_equal(length(testPlt1$layout$facet_params$rows), 1)
    expect_equal(length(testPlt1$layout$facet_params$cols), 1)
    expect_equal(length(testPlt2$layout$facet_params$rows), 0)
    expect_equal(length(testPlt2$layout$facet_params$cols), 1)

})


test_that("plotManhattan utilities work", {
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

    ## Test data 1 -----
    testPrime <- primeManhattanData(
        list(
            "assocStats" = tableReport(fastRep, "FastAssociation"),
            "trait"      = NULL,
            "threshold"  = 2
        )
    )

    ## Check object type -----
    expect_true(is(testPrime, "data.frame"))

    ## Check threshold ----
    threshDf <- testPrime[testPrime$highlight_flag, ]
    expect_false(any(threshDf$p > 0.01))

    ## Check error
    expect_error(
        object = primeManhattanData(
            list(
                "assocStats" = mtcars,
                "trait"      = NULL,
                "threshold"  = 2
            )
        ),
        regexp = "'Chr' column not found in stats dataframe"
    )
    expect_error(
        object = primeManhattanData(
            list(
                "assocStats" = data.frame(
                    "Chr" = c("1", "2", "3"),
                    "Pos" = c(100, 200, 300),
                    "not_important" = c("a", "b", "c")
                ),
                "trait"     = NULL,
                "threshold" = 2
            )
        ),
        regexp = "'Trait' column not found in stats dataframe"
    )
    expect_error(
        object = primeManhattanData(
            list(
                "assocStats" = data.frame(
                    "Chr" = c("1", "2", "3"),
                    "Trait" = c(
                        "important_trait_1",
                        "more_important_trait_2",
                        "the_most_important_trait_in_known_existence"
                    ),
                    "not_important" = c("a", "b", "c")
                ),
                "trait" = NULL,
                "threshold" = 2
            )
        ),
        regexp = "'Pos' column not found in stats dataframe"
    )

})


