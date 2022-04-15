# === Tests for `linkageDiseq()` ====================================

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



## Error tests ----

### General errors ----
test_that("linkageDiseq() throws general exceptions.", {
    expect_error(
        object = linkageDiseq(
            tasObj     = mtcars,
            ldType     = "slidingWindow",
            windowSize = NULL,
            hetCalls   = "missing"
        ),
        regexp = "tasObj is not of class \"TasselGenotypePhenotype\""
    )
    expect_error(
        object = linkageDiseq(
            tasObj     = tasPheno,
            ldType     = "slidingWindow",
            windowSize = NULL,
            hetCalls   = "missing"
        ),
        regexp = "tasObj does contain a Genotype object"
    )
    expect_that(
        object = linkageDiseq(
            tasObj     = tasGenoPhenoFast,
            ldType     = "All",
            windowSize = NULL,
            hetCalls   = "mising"
        ),
        condition = throws_error()
    )
    expect_that(
        object = linkageDiseq(
            tasObj     = tasGenoPhenoFast,
            ldType     = "Everything",
            windowSize = NULL,
            hetCalls   = "missing"
        ),
        condition = throws_error()
    )

})



## Return tests ----
test_that("linkageDiseq() returns correct data structure.", {
    ldDF <- linkageDiseq(
        tasObj     = tasGeno,
        ldType     = "SlidingWindow",
        windowSize = 10,
        hetCalls   = "missing"
    )

    expect_equal(dim(ldDF), c(30875, 17))

    expect_equal(
        object   = colnames(ldDF),
        expected = c(
            "Locus1", "Position1", "Site1", "NumberOfStates1", "States1",
            "Frequency1", "Locus2", "Position2", "Site2", "NumberOfStates2",
            "States2", "Frequency2", "Dist_bp", "R^2", "DPrime", "pDiseq", "N"
        )
    )
})


