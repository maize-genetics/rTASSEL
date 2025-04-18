# === Tests for association plotting utilities ======================

test_that("utility methods work correctly", {
    ### Start logging info
    startLogger()

    ### Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )

    ### Read data - need only non missing data!
    phenoPathFast <- system.file(
        "extdata",
        "mdp_traits_nomissing.txt",
        package = "rTASSEL"
    )

    ### Create rTASSEL phenotype only object
    tasPheno <- readPhenotypeFromPath(
        path = phenoPathFast
    )

    ### Create rTASSEL genotype only object
    tasGeno <- readGenotypeTableFromPath(
        path = genoPathHMP
    )

    ### Create rTASSEL object - use prior TASSEL genotype object
    tasGenoPhenoFast <- readGenotypePhenotype(
        genoPathOrObj = genoPathHMP,
        phenoPathDFOrObj = phenoPathFast
    )

    fastRep <- rTASSEL::assocModelFitter(
        tasGenoPhenoFast,
        . ~ .,
        fitMarkers = TRUE,
        maxP = 1
    )
    tasBLUE <- rTASSEL::assocModelFitter(
        tasGenoPhenoFast,
        . ~ .,
        fitMarkers = FALSE
    )

    expect_warning(
        object = traitValidityChecker(c("dpoll", "Earrrr"), fastRep),
        regexp = "Some traits not found in results and will not be plotted"
    )

})

