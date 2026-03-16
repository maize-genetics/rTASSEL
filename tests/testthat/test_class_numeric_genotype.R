gtSm <- readGenotype(rtMatrices$num_gt_sm)
gtMd <- readGenotype(rtMatrices$num_gt_md)
gtLg <- readGenotype(rtMatrices$num_gt_lg)


# === Basic class structure ========================================

test_that("basic properties", {
    expect_true(is(gtSm, "TasselNumericGenotype"))
    expect_true(is(gtLg, "TasselNumericGenotype"))

    expect_s4_class(gtSm, "TasselGenotype")
    expect_s4_class(gtSm@jRefObj, "jobjRef")
    expect_type(gtSm@jMemAddress, "character")
    expect_type(gtSm@jClass, "character")
    expect_true(is(gtSm@jRefObj, "jobjRef"))
})

test_that("TasselNumericGenotype inherits from TasselGenotype", {
    expect_true(is(gtSm, "TasselGenotype"))
    expect_true(extends("TasselNumericGenotype", "TasselGenotype"))
})


# === show method ==================================================

# Test show method displays correct output
test_that("show method calls with correct parameters", {
    expect_output(show(gtSm), "TasselNumericGenotype")
    expect_output(show(gtLg), "TasselNumericGenotype")
    expect_output(show(gtLg), cli::symbol$ellipsis)
})

test_that("show method prints taxa count, site count, and memory address", {
    out <- capture.output(show(gtSm))
    combined <- paste(out, collapse = "\n")

    expect_match(combined, "taxa")
    expect_match(combined, "sites")
    expect_match(combined, "0x")
})

test_that("show method displays correct dimensions for small matrix", {
    out <- capture.output(show(gtSm))
    combined <- paste(out, collapse = "\n")

    nTaxa  <- gtSm@jRefObj$numberOfTaxa()
    nSites <- gtSm@jRefObj$numberOfSites()

    expect_match(combined, as.character(nTaxa),  fixed = TRUE)
    expect_match(combined, as.character(nSites), fixed = TRUE)
})


# === Slot content validation ======================================

test_that("jMemAddress slot is a hex-like string", {
    expect_match(gtSm@jMemAddress, "^[0-9a-f]+$", ignore.case = TRUE)
    expect_match(gtLg@jMemAddress, "^[0-9a-f]+$", ignore.case = TRUE)
})

test_that("jClass slot contains a Java class path", {
    expect_match(gtSm@jClass, "\\.", perl = TRUE)
})

test_that("dispData slot is a non-empty list", {
    expect_true(length(gtSm@dispData) > 0)
    expect_true(length(gtLg@dispData) > 0)
})

test_that("different matrix sizes produce distinct memory addresses", {
    expect_false(identical(gtSm@jMemAddress, gtLg@jMemAddress))
})


# === javaRefObj method ============================================

test_that("javaRefObj returns the exact jRefObj slot", {
    jRef <- javaRefObj(gtSm)
    expect_s4_class(jRef, "jobjRef")
    expect_identical(jRef, gtSm@jRefObj)
})

test_that("javaRefObj result has callable methods", {
    jRef <- javaRefObj(gtSm)
    expect_true(jRef$numberOfTaxa() > 0)
    expect_true(jRef$numberOfSites() > 0)
})

test_that("javaRefObj taxa/site counts match input matrix", {
    jRef <- javaRefObj(gtMd)
    expect_equal(jRef$numberOfTaxa(),  nrow(rtMatrices$num_gt_md))
    expect_equal(jRef$numberOfSites(), ncol(rtMatrices$num_gt_md))
})


# === Construction from varying matrix sizes =======================

test_that("readGenotype produces correct taxa count from matrix", {
    expect_equal(gtSm@jRefObj$numberOfTaxa(), nrow(rtMatrices$num_gt_sm))
    expect_equal(gtMd@jRefObj$numberOfTaxa(), nrow(rtMatrices$num_gt_md))
    expect_equal(gtLg@jRefObj$numberOfTaxa(), nrow(rtMatrices$num_gt_lg))
})

test_that("readGenotype produces correct site count from matrix", {
    expect_equal(gtSm@jRefObj$numberOfSites(), ncol(rtMatrices$num_gt_sm))
    expect_equal(gtMd@jRefObj$numberOfSites(), ncol(rtMatrices$num_gt_md))
    expect_equal(gtLg@jRefObj$numberOfSites(), ncol(rtMatrices$num_gt_lg))
})


# === Matrix validation errors =====================================

test_that("readNumericGenotypeFromRMatrix rejects non-matrix input", {
    expect_error(
        readNumericGenotypeFromRMatrix(data.frame(x = 1)),
        "not of type 'matrix'"
    )
})

test_that("readNumericGenotypeFromRMatrix rejects matrix without row names", {
    m <- matrix(1:4, nrow = 2)
    colnames(m) <- c("s1", "s2")
    expect_error(
        readNumericGenotypeFromRMatrix(m),
        "No row names"
    )
})

test_that("readNumericGenotypeFromRMatrix rejects matrix without column names", {
    m <- matrix(1:4, nrow = 2)
    rownames(m) <- c("t1", "t2")
    expect_error(
        readNumericGenotypeFromRMatrix(m),
        "No column names"
    )
})


