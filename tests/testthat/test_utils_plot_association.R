# === Tests for association plotting utilities ======================

test_that("utility methods work correctly", {
    ### Start logging info
    startLogger()

    ### Shared fixtures (see helper_vars.R)
    tasDataset <- rtObjs$ds_hmp_ph_nomiss

    fastRep <- rTASSEL::assocModelFitter(
        tasDataset,
        . ~ .,
        fitMarkers = TRUE,
        maxP = 1
    )
    tasBLUE <- rTASSEL::assocModelFitter(
        tasDataset,
        . ~ .,
        fitMarkers = FALSE
    )

    expect_warning(
        object = traitValidityChecker(c("dpoll", "Earrrr"), fastRep),
        regexp = "Some traits not found in results and will not be plotted"
    )

})

