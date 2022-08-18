# === Tests SummarizedExperiment creation method ====================

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
tasGenoPheno <- readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoPathFast
)


## Error tests ----
test_that("getSumExpFromGenotypeTable() throws general exceptions.", {
    expect_error(
        object = getSumExpFromGenotypeTable(
            tasObj            = mtcars,
            coerceDosageToInt = FALSE,
            verbose           = FALSE
        ),
        regexp = "`tasObj` must be of class `TasselGenotypePhenotype`"
    )
    expect_error(
        object = getSumExpFromGenotypeTable(
            tasObj            = tasPheno,
            coerceDosageToInt = FALSE,
            verbose           = FALSE
        ),
        regexp = "TASSEL genotype object not found"
    )
})


## Return tests ----
test_that("getSumExpFromGenotypeTable() returns correct data.", {
    tasSE <- getSumExpFromGenotypeTable(
        tasObj            = tasGenoPheno,
        coerceDosageToInt = FALSE,
        verbose           = FALSE
    )

    expect_equal(
        object = class(tasSE)[1],
        expected = "RangedSummarizedExperiment"
    )
    expect_equal(
        object = dim(tasSE),
        expected = c(3093, 278)
    )
})

test_that("getSumExpFromGenotypeTable() returns correct dosage types.", {
    tasSERaw <- getSumExpFromGenotypeTable(
        tasObj            = tasGenoPheno,
        coerceDosageToInt = FALSE,
        verbose           = FALSE
    )
    assaySERaw <- SummarizedExperiment::assay(tasSERaw)[[1, 1]]

    tasSEInt <- getSumExpFromGenotypeTable(
        tasObj            = tasGenoPheno,
        coerceDosageToInt = TRUE,
        verbose           = FALSE
    )
    assaySEInt <- SummarizedExperiment::assay(tasSEInt)[[1, 1]]

    expect_equal(
        object = class(assaySEInt),
        expected = "integer"
    )
    expect_equal(
        object = class(assaySERaw),
        expected = "raw"
    )
})


