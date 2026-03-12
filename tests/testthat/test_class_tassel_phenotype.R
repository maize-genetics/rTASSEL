# === Shared test fixtures ==========================================

phTestDf <- data.frame(
    taxa_id      = c("line_a", "line_b"),
    plant_height = c(12.3, 22.8),
    PC1          = c(0.5, -1.5),
    yield        = c(2, 3),
    stringsAsFactors = FALSE
)
phTestAttr <- data.frame(
    col_id      = c("taxa_id", "plant_height", "PC1", "yield"),
    tassel_attr = c("taxa", "data", "covariate", "data"),
    stringsAsFactors = FALSE
)
phTestObj <- readPhenotype(phTestDf, attr = phTestAttr)


# === TasselPhenotype construction =================================

test_that("TasselPhenotype creation works", {
    expect_no_error(phenotype <- readPhenotype(phTestDf, attr = phTestAttr))
    expect_s4_class(phenotype, "TasselPhenotype")
    expect_true(all(
        c("attrData", "attrSummary", "dispData", "rData",
          "jRefObj", "jMemAddress", "jClass") %in% slotNames(phenotype)
    ))
})


# === readPhenotype dispatch =======================================

test_that("TasselPhenotype methods work correctly", {
    expect_s3_class(attributeData(phTestObj), "data.frame")
    expect_equal(nrow(attributeData(phTestObj)), 4)
    expect_equal(ncol(attributeData(phTestObj)), 5)
    expect_true(inherits(javaRefObj(phTestObj), "jobjRef"))
})

test_that("TasselPhenotype validation works", {
    expect_error(
        readPhenotype(data.frame(x = 1), attr = NULL),
        "Phenotype objects evaluated from 'data.frame' need attribute metadata"
    )
    expect_error(
        readPhenotype(list(x = 1)),
        "Unsupported input type for 'x'"
    )
})

test_that("readPhenotype rejects unsupported types", {
    expect_error(readPhenotype(42),   "Unsupported input type")
    expect_error(readPhenotype(TRUE), "Unsupported input type")
    expect_error(readPhenotype(NULL), "Unsupported input type")
})

test_that("readPhenotype dispatches to file path when given a character", {
    phPath <- readPhenotype(rtFiles$ph_nomiss_path)
    expect_s4_class(phPath, "TasselPhenotype")
    expect_true(nrow(phPath@rData) > 0)
})

test_that("readPhenotype dispatches to data.frame path with attr", {
    ph <- readPhenotype(phTestDf, attr = phTestAttr)
    expect_s4_class(ph, "TasselPhenotype")
    expect_true("plant_height" %in% attributeData(ph)$trait_id)
})


# === Display methods ==============================================

test_that("TasselPhenotype display methods work", {
    expect_output(show(phTestObj), "TasselPhenotype")
    expect_true(all(
        c("nTaxa", "nTraits", "nDfRow", "nCap", "jMem") %in%
            names(attributes(phTestObj@dispData))
    ))
})

test_that("show method prints taxa and trait counts", {
    out <- capture.output(show(phTestObj))
    combined <- paste(out, collapse = "\n")

    expect_match(combined, "taxa")
    expect_match(combined, "traits")
})

test_that("tbl_format_header returns styled header string", {
    hdr <- pillar::tbl_format_header(phTestObj@dispData)
    hdrText <- paste(hdr, collapse = " ")

    expect_match(hdrText, "TasselPhenotype", fixed = TRUE)
    expect_match(hdrText, "taxa")
    expect_match(hdrText, "traits")
})

test_that("printed output includes Java memory address", {
    out <- capture.output(print(phTestObj@dispData))
    combined <- paste(out, collapse = "\n")

    expect_match(combined, "0x", fixed = TRUE)
})


# === Slot content validation ======================================

test_that("attrData contains expected columns", {
    ad <- phTestObj@attrData

    expect_s3_class(ad, "data.frame")
    expect_true("trait_id" %in% colnames(ad))
    expect_true("trait_type" %in% colnames(ad))
    expect_true("taxa" %in% ad$trait_type)
})

test_that("attrSummary is a named list reflecting trait type counts", {
    as <- phTestObj@attrSummary

    expect_type(as, "list")
    expect_true(length(as) > 0)
    expect_true("taxa" %in% names(as))
    expect_equal(as[["taxa"]], 1L)
})

test_that("rData slot is a data.frame with correct dimensions", {
    rd <- phTestObj@rData

    expect_s3_class(rd, "data.frame")
    expect_equal(nrow(rd), nrow(phTestDf))
    expect_true(ncol(rd) >= ncol(phTestDf))
})

test_that("jMemAddress slot is a hex-like string", {
    expect_match(phTestObj@jMemAddress, "^[0-9a-f]+$", ignore.case = TRUE)
})

test_that("jClass slot contains a Java class path", {
    expect_match(phTestObj@jClass, "\\.", perl = TRUE)
})

test_that("dispData slot is a java_pheno_tbl", {
    expect_s3_class(phTestObj@dispData, "java_pheno_tbl")
})


# === javaRefObj method ============================================

test_that("javaRefObj returns the exact jRefObj slot", {
    jRef <- javaRefObj(phTestObj)
    expect_s4_class(jRef, "jobjRef")
    expect_identical(jRef, phTestObj@jRefObj)
})

test_that("javaRefObj result exposes taxa count", {
    jRef <- javaRefObj(phTestObj)
    expect_true(jRef$taxa()$numberOfTaxa() > 0)
})


# === attributeData method =========================================

test_that("attributeData returns same object as attrData slot", {
    expect_identical(attributeData(phTestObj), phTestObj@attrData)
})

test_that("attributeData trait_id values match input column names", {
    traitIds <- attributeData(phTestObj)$trait_id
    expect_true(all(c("plant_height", "PC1", "yield") %in% traitIds))
})


# === pillar_shaft S3 methods ======================================

test_that("pillar_shaft.taxa produces a pillar shaft", {
    x <- taxaVctr(c("line_a", "line_b"))
    shaft <- pillar::pillar_shaft(x)
    expect_s3_class(shaft, "pillar_shaft")
})

test_that("pillar_shaft.data produces a pillar shaft", {
    x <- dataVctr(c(1.1, 2.2))
    shaft <- pillar::pillar_shaft(x)
    expect_s3_class(shaft, "pillar_shaft")
})

test_that("pillar_shaft.cov produces a pillar shaft", {
    x <- covVctr(c(0.5, -1.5))
    shaft <- pillar::pillar_shaft(x)
    expect_s3_class(shaft, "pillar_shaft")
})

test_that("pillar_shaft.fact produces a pillar shaft", {
    x <- factVctr(c("A", "B"))
    shaft <- pillar::pillar_shaft(x)
    expect_s3_class(shaft, "pillar_shaft")
})


# === File-based phenotype construction ============================

test_that("readPhenotype from file produces consistent slots", {
    ph <- readPhenotype(rtFiles$ph_nomiss_path)

    expect_s4_class(ph, "TasselPhenotype")
    expect_s3_class(ph@rData, "data.frame")
    expect_s3_class(ph@attrData, "data.frame")
    expect_type(ph@attrSummary, "list")
    expect_true(nzchar(ph@jMemAddress))
    expect_true(nzchar(ph@jClass))
})

test_that("readPhenotype from file matches readPhenotypeFromPath trait names", {
    phNew <- readPhenotype(rtFiles$ph_nomiss_path)
    phOld <- rtObjs$ph_nomiss

    expect_equal(
        names(getPhenotypeDF(phOld)),
        colnames(phNew@rData)
    )
})


