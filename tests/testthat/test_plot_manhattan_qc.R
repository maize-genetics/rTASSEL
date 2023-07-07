# === Tests for Manhattan QC plotting ===============================

test_that("plotManhattanQC works correctly", {
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
        fitMarkers = TRUE,
        maxP = 1
    )

    testGRData <- GenomicRanges::GRanges(
        seqnames = c("8", "2", "5"),
        ranges = IRanges::IRanges(
            start = c(131576889, 3000000, 1000000),
            end   = c(131580316, 3070000, 1050000)
        ),
        range_id = c("Zm01", "Zm02", "Zm03"),
        classical_id = c("rap2.7", "classic_02", "classic_03")
    )
    testPlt <- plotManhattanQC(fastRep, gr = testGRData)
    expect_true(is(testPlt, "gg"))
    expect_error(
        object = plotManhattanQC(fastRep, gr = testGRData, window = 100),
        regexp = "No association results found for any trait"
    )

    testGRData <- GenomicRanges::GRanges(
        seqnames = c("8", "2", "5"),
        ranges = IRanges::IRanges(
            start = c(131576889, 3000000, 1000000),
            end   = c(131580316, 3070000, 1050000)
        ),
        range_id = c("Zm01", "Zm02", "Zm03")
    )
    expect_warning(
        object = plotManhattanQC(fastRep, gr = testGRData, classicNames = TRUE),
        regexp = "'classical_id' parameter not found."
    )
    expect_message(
        object = plotManhattanQC(fastRep, gr = testGRData),
        regexp = "No association results found"
    )


    tasBLUE <- rTASSEL::assocModelFitter(
        tasGenoPhenoFast,
        . ~ .,
        fitMarkers = FALSE
    )
    expect_error(plotManhattanQC(mtcars, gr = testGRData))
    expect_error(plotManhattanQC(tasBLUE, gr = testGRData))
})

