# === Tests for genotype table functions ============================

## Preamble - load data ----

### Start logging info
startLogger()

### Load hapmap data
genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)

### Read data - need only non missing data!
phenoPathFast <- system.file(
    "extdata",
    "mdp_traits_nomissing.txt",
    package = "rTASSEL"
)

### Create rTASSEL phenotype only object
tasPheno <- readPhenotypeFromPath(
    path = phenoPathFast
)

### Create rTASSEL genotype only object
tasGeno <- readGenotypeTableFromPath(
    path = genoPathHMP
)

### Create rTASSEL object - use prior TASSEL genotype object
tasGenoPhenoFast <- readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoPathFast
)

### Filter object for further tests
filterGenoObj <- filterGenotypeTableSites(
    tasObj = tasGeno,
    siteRangeFilterType = "sites",
    startSite = 0,
    endSite = 10
)
filterGenoObj <- filterGenotypeTableTaxa(
    tasObj = filterGenoObj,
    taxa = taxaList(tasGeno)[grep("^[0-9]|^A", taxaList(tasGeno))]
)


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
    expect_error(as.matrix.TasselGenotypePhenotype(tasPheno))
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


