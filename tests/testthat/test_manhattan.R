# === Tests for Manhattan Plot functions ============================

## Error tests ----
test_that("Manhattan Plot throws error.", {

    ## Load data ----

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

    ### Create rTASSEL object - use prior TASSEL genotype object
    tasGenoPhenoFast <- rTASSEL::readGenotypePhenotype(
        genoPathOrObj = genoPathHMP,
        phenoPathDFOrObj = phenoPathFast
    )


    ## Fast association ----
    fastRep <- rTASSEL::assocModelFitter(
        tasGenoPhenoFast,
        . ~ .,
        fastAssociation = TRUE,
        fitMarkers = TRUE,
        maxP = 1
    )

    expect_error(
        manhattanPlot(
            assocStats = fastRep$FastAssociation,
            trait = "dpoll",
            showSigMarkers = TRUE,
            showRug = TRUE,
            colors = c("lightgrey", "red", "blue")
        ),
        "Please enter a numeric threshold value."
    )

})


## Equality tests ----
test_that("Manhattan Plot returns correct plot layers", {

    ## Load data ----

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

    ### Create rTASSEL object - use prior TASSEL genotype object
    tasGenoPhenoFast <- rTASSEL::readGenotypePhenotype(
        genoPathOrObj = genoPathHMP,
        phenoPathDFOrObj = phenoPathFast
    )


    ## Fast association ----
    fastRep <- rTASSEL::assocModelFitter(
        tasGenoPhenoFast,
        . ~ .,
        fastAssociation = TRUE,
        fitMarkers = TRUE,
        maxP = 1
    )

    ## Plot object ----
    p <- manhattanPlot(
        assocStats     = fastRep$FastAssociation,
        trait          = "EarHT",
        showSigMarkers = TRUE,
        showRug        = TRUE,
        threshold      = 6
    )


    ## Tests ----
    expect_s3_class(
        object = p$layers[[1]],
        class  = "ggproto"
    )
    expect_s3_class(
        object = p$layers[[1]]$geom,
        class  = "GeomPoint"
    )
    expect_s3_class(
        object = p$layers[[3]]$geom,
        class  = "GeomRug"
    )
    expect_equal(
        object   = p$labels$title,
        expected = "Trait: EarHT"
    )
    expect_equal(
        object   = p$labels$x,
        expected = "Position"
    )
    expect_equal(
        object   = p$labels$y,
        expected = ~-log[10] ~ "(" * italic(p) * "-value)"
    )
})


