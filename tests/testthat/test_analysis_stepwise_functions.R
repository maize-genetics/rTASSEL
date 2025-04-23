test_that("stepwiseModelFitter basic operations", {
    expect_error(stepwiseModelFitter(mtcars))
    expect_error(stepwiseModelFitter(rtObjs$gt_hmp_ph_nomiss, entryLimit = -0.1))
    expect_error(stepwiseModelFitter(rtObjs$gt_hmp_ph_nomiss, exitLimit = 1.5))
    expect_error(stepwiseModelFitter(rtObjs$gt_hmp_ph_nomiss, maxNumberOfMarkers = 1e5))
    expect_error(stepwiseModelFitter(rtObjs$gt_hmp_ph_nomiss, nPermutations = 1e6))

    gtSmall <- filterGenotypeTableSites(
        rtObjs$gt_hmp,
        siteRangeFilterType = "sites",
        startSite = 0, endSite = 100
    )
    gtPhSmall <- readGenotypePhenotype(gtSmall, rtObjs$ph_nomiss)

    stepRes01 <- stepwiseModelFitter(gtPhSmall)
    stepRes02 <- stepwiseModelFitter(gtPhSmall, dpoll ~ .)

    expect_true(is(stepRes01, "AssociationResults"))
    expect_equal(traitNames(stepRes02), "dpoll")

    # Marker validation
    # NOTE: marker set validated with:
    #         * TASSEL 5 GUI
    #         * external OLS libraries in R ("olsrr")
    stepRes02TabRep <- tableReport(stepRes02)
    truthSet <- c(
        "PZA00731.6", "PZD00098.1", "PHM4531.46",
        "PZA00447.8", "PZA00181.2", "PZA02487.1",
        "PZA03128.3", "PZA02921.4", "PZA00258.3"
    )
    expect_true(all(truthSet %in% stepRes02TabRep$Name))
})


