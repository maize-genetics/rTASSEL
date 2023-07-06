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
        seqnames = c("1", "2", "5"),
        ranges = IRanges::IRanges(
            start = c(1,      3000000, 1000000),
            end   = c(250000, 3005000, 1050000)
        ),
        range_id = c("Zm01", "Zm02", "Zm03"),
        canonical_id = c("Canon_01", "Canon_02", "Canon_03")
    )

    testPlt <- plotManhattanQC(fastRep, testGRData, window = 1000000)
    expect_true(is(testPlt, "gg"))
})

