# === Tests for genotype table functions ============================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasPheno   <- rtObjs$ph_nomiss
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss

### A small corner of chromosome 1 for the summary tests
filterGenoObj <- tasGeno[
    taxaWhere(grepl("^[0-9]|^A", taxaId)),
    sites(1:11)
]


test_that("readGenotypeTableFromPath()", {
    expect_error(readGenotypeTableFromPath("i/dont/exist"))
})

test_that("getSumExpFromGenotypeTable()", {
    expect_message(getSumExpFromGenotypeTable(filterGenoObj, verbose = TRUE))
})

test_that("getMinMaxPhysPositions()", {
    expect_error(getMinMaxPhysPositions(mtcars))
    expect_error(getMinMaxPhysPositions(tasPheno))

    testObj <- getMinMaxPhysPositions(filterGenoObj)
    expect_equal(names(testObj), "1")
    expect_equal(length(testObj[[1]]), 2)

    testObj <- getMinMaxPhysPositions(tasGeno)
    expect_gt(length(testObj), 1)
    expect_equal(length(testObj[[1]]), 2)
})

test_that("getMinMaxVarSites()", {
    expect_error(getMinMaxVarSites(mtcars))
    expect_error(getMinMaxVarSites(tasPheno))

    testObj <- getMinMaxVarSites(filterGenoObj)
    expect_equal(names(testObj), "1")
    expect_equal(length(testObj[[1]]), 2)

    testObj <- getMinMaxVarSites(tasGeno)
    expect_gt(length(testObj), 1)
    expect_equal(length(testObj[[1]]), 2)
})

test_that("as.matrix.TasselGenotypePhenotyp()", {
    expect_error(as.matrix.TasselGenotypePhenotype(mtcars))
    expect_error(as.matrix.TasselGenotypePhenotype(rtObjsLegacy$ph_nomiss))
})

test_that("siteSummary()", {
    expect_error(siteSummary(mtcars))
    expect_error(siteSummary(tasPheno))

    testObj <- siteSummary(filterGenoObj)
    truthColNames <- c(
        "Site_Number",
        "Site_Name",
        "Chromosome",
        "Physical_Position",
        "Number_of_Taxa",
        "Ref",
        "Alt",
        "Major_Allele",
        "Major_Allele_Gametes",
        "Major_Allele_Proportion",
        "Major_Allele_Frequency",
        "Minor_Allele",
        "Minor_Allele_Gametes",
        "Minor_Allele_Proportion",
        "Minor_Allele_Frequency",
        "Allele_3",
        "Allele_3_Gametes",
        "Allele_3_Proportion",
        "Allele_3_Frequency",
        "Allele_4",
        "Allele_4_Gametes",
        "Allele_4_Proportion",
        "Allele_4_Frequency",
        "Allele_5",
        "Allele_5_Gametes",
        "Allele_5_Proportion",
        "Allele_5_Frequency",
        "Allele_6",
        "Allele_6_Gametes",
        "Allele_6_Proportion",
        "Allele_6_Frequency",
        "Allele_7",
        "Allele_7_Gametes",
        "Allele_7_Proportion",
        "Allele_7_Frequency",
        "Gametes_Missing",
        "Proportion_Missing",
        "Number_Heterozygous",
        "Proportion_Heterozygous",
        "Inbreeding_Coefficient",
        "Inbreeding_Coefficient_Scaled_by_Missing"
    )
    expect_equal(nrow(testObj), 11)
    expect_equivalent(colnames(testObj), truthColNames)
    expect_true(inherits(siteSummary(filterGenoObj), "data.frame"))
})

test_that("taxaSummary()", {
    expect_error(taxaSummary(mtcars))
    expect_error(taxaSummary(tasPheno))

    testObj <- taxaSummary(filterGenoObj)
    truthColNames <- c(
        "Taxa",
        "Taxa_Name",
        "Number_of_Sites",
        "Gametes_Missing",
        "Proportion_Missing",
        "Number_Heterozygous",
        "Proportion_Heterozygous",
        "Inbreeding_Coefficient",
        "Inbreeding_Coefficient_Scaled_by_Missing"
    )
    expect_equal(nrow(testObj), 24)
    expect_equivalent(colnames(testObj), truthColNames)
    expect_true(inherits(siteSummary(filterGenoObj), "data.frame"))
})

test_that("summaries and position ranges accept a genomic dataset", {
    # A dataset summarises its joined genotype table, which has dropped the
    # taxa without phenotype records
    expect_equal(siteSummary(tasDataset), siteSummary(genotype(tasDataset)))
    expect_equal(nrow(taxaSummary(tasDataset)), 278)
    expect_equal(
        getMinMaxPhysPositions(tasDataset),
        getMinMaxPhysPositions(tasGeno)
    )
    expect_equal(getMinMaxVarSites(tasDataset), getMinMaxVarSites(tasGeno))
})


## Back-compatibility ----
test_that("genotype table functions accept a deprecated TasselGenotypePhenotype", {
    legacyGt <- rtObjsLegacy$gt_hmp

    expect_equal(siteSummary(legacyGt), siteSummary(tasGeno))
    expect_equal(taxaSummary(legacyGt), taxaSummary(tasGeno))
    expect_equal(
        getMinMaxPhysPositions(legacyGt),
        getMinMaxPhysPositions(tasGeno)
    )
    expect_equal(getMinMaxVarSites(legacyGt), getMinMaxVarSites(tasGeno))
})


