# === Tests for association results sub-classes =====================

test_that("Method dispatch for sub classes works.", {
    ## Preamble - load data ----
    startLogger()

    tasPheno   <- rtObjs$ph_nomiss
    tasDataset <- rtObjs$ds_hmp_ph_nomiss
    tasKin     <- kinshipMatrix(tasDataset)

    ### Association objects
    tasBLUE <- assocModelFitter(
        tasObj  = tasPheno,
        formula = . ~ .
    )
    tasGLM <- assocModelFitter(
        tasObj     = tasDataset,
        formula    = . ~ .,
        fitMarkers = TRUE
    )
    tasMLM <- assocModelFitter(
        tasObj     = tasDataset,
        formula    = . ~ .,
        fitMarkers = TRUE,
        kinship    = tasKin
    )
    tasFast <- assocModelFitter(
        tasObj          = tasDataset,
        formula         = . ~ .,
        fitMarkers      = TRUE,
        fastAssociation = TRUE
    )

    ## Tests ----
    dimBLUE <- dim(tableReport(tasBLUE, "BLUE"))
    dimGLM  <- dim(tableReport(tasGLM, "GLM_Stats"))
    dimMLM  <- dim(tableReport(tasMLM, "MLM_Stats"))
    dimFast <- dim(tableReport(tasFast, "FastAssociation"))

    expect_true(is(tableReport(tasBLUE), "data.frame"))
    expect_true(is(tableReport(tasGLM), "data.frame"))
    expect_true(is(tableReport(tasMLM), "data.frame"))
    expect_true(is(tableReport(tasFast), "data.frame"))

    expect_equal(dim(tableReport(tasBLUE)), dimBLUE)
    expect_equal(dim(tableReport(tasGLM)), dimGLM)
    expect_equal(dim(tableReport(tasMLM)), dimMLM)
    expect_equal(dim(tableReport(tasFast)), dimFast)
})
