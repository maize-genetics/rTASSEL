# === Tests for loading demo data ===================================

test_that("loadDemoData loads phenotype data correctly", {
    phenoData <- loadDemoData("phenotype")
    expect_true(!is.null(phenoData))
    expect_true(is(phenoData, "TasselGenotypePhenotype"))
})

test_that("loadDemoData loads genotype data correctly", {
    genoData <- loadDemoData("genotype")
    expect_true(!is.null(genoData))
    expect_true(inherits(genoData, "TasselGenotypePhenotype"))
})

test_that("loadDemoData loads combined data correctly", {
    combinedData <- loadDemoData("combine")
    expect_true(!is.null(combinedData))
    expect_true(inherits(combinedData, "TasselGenotypePhenotype"))
})

test_that("loadDemoData handles invalid type correctly", {
    expect_error(loadDemoData("invalid"))
})


