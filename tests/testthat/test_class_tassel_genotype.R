
test_that("TasselGenotype class construction and methods work", {

    gtTest <- readGenotype(rtFiles$gt_hmp_path)


    # Test class structure
    expect_s4_class(gtTest, "TasselGenotype")
    expect_type(gtTest@dispData, "list")
    expect_s4_class(gtTest@jRefObj, "jobjRef")
    expect_type(gtTest@jMemAddress, "character")
    expect_type(gtTest@jClass, "character")
    expect_true(is(gtTest@jRefObj, "jobjRef"))

    # Test show method (capturing output)
    expect_output(show(gtTest), "TasselGenotype")  # Basic check that output contains class name
})

test_that("readGenotype handles invalid inputs correctly", {
    # Test invalid file path
    expect_error(readGenotype("nonexistent/path.txt"))

    # Test invalid input type
    expect_error(readGenotype(list()), "Unsupported data type")
})
