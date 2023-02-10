# === Tests for taxa lists ==========================================

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

## Tests ----
test_that("getTaxaList() returns correct data", {
    testObj <- getTaxaList(mtcars)
    expect_true(rJava::is.jnull(testObj))

    testObj <- getTaxaList(tasGenoPhenoFast@jPositionList)
    expect_true(rJava::is.jnull(testObj))
})

test_that("taxaList() returns correct excpetions", {
    expect_error(taxaList(mtcars))
})

test_that("getTaxaIDs() returns correct excpetions", {
    expect_error(getTaxaIDs(mtcars))
})



