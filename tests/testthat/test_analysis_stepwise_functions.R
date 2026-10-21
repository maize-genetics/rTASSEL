test_that("stepwiseModelFitter basic operations", {
    expect_error(stepwiseModelFitter(mtcars))
    expect_error(stepwiseModelFitter(rtObjs$ds_hmp_ph_nomiss, entryLimit = -0.1))
    expect_error(stepwiseModelFitter(rtObjs$ds_hmp_ph_nomiss, exitLimit = 1.5))
    expect_error(stepwiseModelFitter(rtObjs$ds_hmp_ph_nomiss, maxNumberOfMarkers = 1e5))
    expect_error(stepwiseModelFitter(rtObjs$ds_hmp_ph_nomiss, nPermutations = 1e6))

    dsSmall <- rtObjs$ds_hmp_ph_nomiss[, sites(1:101)]

    stepRes01 <- stepwiseModelFitter(dsSmall)
    stepRes02 <- stepwiseModelFitter(dsSmall, dpoll ~ .)

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


test_that("stepwiseModelFitter accepts a deprecated TasselGenotypePhenotype", {
    legacySmall <- filterGenotypeTableSites(
        rtObjsLegacy$gt_hmp_ph_nomiss,
        siteRangeFilterType = "sites",
        startSite = 0, endSite = 100
    )

    expect_equal(
        tableReport(stepwiseModelFitter(legacySmall, dpoll ~ .)),
        tableReport(
            stepwiseModelFitter(rtObjs$ds_hmp_ph_nomiss[, sites(1:101)], dpoll ~ .)
        )
    )
})
