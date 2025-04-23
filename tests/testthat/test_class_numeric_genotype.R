gtSm <- readGenotype(rtMatrices$num_gt_sm)
gtLg <- readGenotype(rtMatrices$num_gt_lg)

# Test basic characteristics
test_that("basic properties", {
    expect_true(is(gtSm, "TasselNumericGenotype"))
    expect_true(is(gtLg, "TasselNumericGenotype"))

    # Test class structure
    expect_s4_class(gtSm, "TasselGenotype")
    expect_type(gtSm@dispData, "list")
    expect_s4_class(gtSm@jRefObj, "jobjRef")
    expect_type(gtSm@jMemAddress, "character")
    expect_type(gtSm@jClass, "character")
    expect_true(is(gtSm@jRefObj, "jobjRef"))
})


# Test show method calls printNumGtDisp with correct parameters
test_that("show method calls with correct parameters", {
    expect_output(show(gtSm), "TasselNumericGenotype")
    expect_output(show(gtLg), "TasselNumericGenotype")
    expect_output(show(gtLg), cli::symbol$ellipsis)
})


