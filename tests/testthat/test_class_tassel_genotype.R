
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


## readGenotype edge cases for unsupported types ----

test_that("readGenotype rejects non-character non-matrix types", {
    expect_error(readGenotype(42),          "Unsupported data type")
    expect_error(readGenotype(TRUE),        "Unsupported data type")
    expect_error(readGenotype(data.frame()), "Unsupported data type")
    expect_error(readGenotype(NULL),        "Unsupported data type")
})

test_that("readGenotype reports correct error for missing file", {
    expect_error(
        readGenotype("/tmp/definitely_not_a_real_genotype_file.hmp.txt"),
        "not a valid file"
    )
})


## readGenotype with matrix input ----

test_that("readGenotype accepts a matrix and returns a genotype object", {
    m <- rtMatrices$num_gt_md
    gt <- readGenotype(m)

    expect_s4_class(gt, "TasselNumericGenotype")
    expect_type(gt@dispData, "list")
    expect_s4_class(gt@jRefObj, "jobjRef")
    expect_true(nzchar(gt@jMemAddress))
})


## readGenotype with sortPositions / keepDepth flags ----

test_that("readGenotype with sortPositions = TRUE produces valid object", {
    gt <- readGenotype(rtFiles$gt_hmp_path, sortPositions = TRUE)
    expect_s4_class(gt, "TasselGenotype")
    expect_true(gt@jRefObj$numberOfTaxa() > 0)
    expect_true(gt@jRefObj$numberOfSites() > 0)
})

test_that("readGenotype with keepDepth = TRUE produces valid object", {
    gt <- readGenotype(rtFiles$gt_vcf_path, keepDepth = TRUE)
    expect_s4_class(gt, "TasselGenotype")
    expect_true(gt@jRefObj$numberOfTaxa() > 0)
})


## show method detail ----

test_that("show method prints taxa count, site count, and memory address", {
    gtTest <- readGenotype(rtFiles$gt_hmp_path)

    out <- capture.output(show(gtTest))
    combined <- paste(out, collapse = "\n")

    expect_match(combined, "taxa")
    expect_match(combined, "sites")
    expect_match(combined, "0x")
})


## javaRefObj method ----

test_that("javaRefObj returns the Java reference for TasselGenotype", {
    gtTest <- readGenotype(rtFiles$gt_hmp_path)

    jRef <- javaRefObj(gtTest)
    expect_s4_class(jRef, "jobjRef")
    expect_identical(jRef, gtTest@jRefObj)
})

test_that("javaRefObj result has callable methods", {
    gtTest <- readGenotype(rtFiles$gt_hmp_path)
    jRef <- javaRefObj(gtTest)

    expect_true(jRef$numberOfTaxa() > 0)
    expect_true(jRef$numberOfSites() > 0)
})


## Slot integrity ----

test_that("jMemAddress slot is a hex-like string", {
    gtTest <- readGenotype(rtFiles$gt_hmp_path)
    expect_match(gtTest@jMemAddress, "^[0-9a-f]+$", ignore.case = TRUE)
})

test_that("jClass slot contains a Java class path", {
    gtTest <- readGenotype(rtFiles$gt_hmp_path)
    expect_match(gtTest@jClass, "\\.", perl = TRUE)
})

test_that("dispData slot is a non-empty list", {
    gtTest <- readGenotype(rtFiles$gt_hmp_path)
    expect_true(length(gtTest@dispData) > 0)
})
