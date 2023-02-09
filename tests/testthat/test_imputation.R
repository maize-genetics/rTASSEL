# === Tests for imputation methods ==================================

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



## Imputation (Numeric) ----
test_that("imputeNumeric() returns correct exceptions", {
    expect_error(
        object = imputeNumeric(mtcars),
        regexp = "`tasObj` must be of class `TasselGenotypePhenotype`"
    )

    expect_error(
        object = imputeNumeric(tasPheno),
        regexp = "TASSEL genotype object not found"
    )
})


test_that("imputeNumeric() returns correct data", {
    expect_true(
        inherits(
            imputeNumeric(
                tasObj = filterGenoObj,
                byMean = TRUE,
                nearestNeighbors = 5,
                distance = "Euclidean"
            ),
            "TasselGenotypePhenotype"
        )
    )
    expect_true(
        inherits(
            imputeNumeric(
                tasObj = filterGenoObj,
                byMean = FALSE,
                nearestNeighbors = 5,
                distance = "Euclidean"
            ),
            "TasselGenotypePhenotype"
        )
    )
    expect_true(
        inherits(
            imputeNumeric(
                tasObj = filterGenoObj,
                byMean = FALSE,
                nearestNeighbors = 5,
                distance = "Manhattan"
            ),
            "TasselGenotypePhenotype"
        )
    )
    expect_true(
        inherits(
            imputeNumeric(
                tasObj = filterGenoObj,
                byMean = FALSE,
                nearestNeighbors = 5,
                distance = "Cosine"
            ),
            "TasselGenotypePhenotype"
        )
    )
})



## Imputation (LD KNNi) ----
test_that("imputeLDKNNi() returns correct exceptions", {
    expect_error(
        object = imputeLDKNNi(mtcars),
        regexp = "`tasObj` must be of class `TasselGenotypePhenotype`"
    )

    expect_error(
        object = imputeLDKNNi(tasPheno),
        regexp = "TASSEL genotype object not found"
    )

    expect_error(object = imputeLDKNNi(filterGenoObj, highLDSSites = 1))

    expect_error(object = imputeLDKNNi(filterGenoObj, highLDSSites = -1))

    expect_error(object = imputeLDKNNi(filterGenoObj, highLDSSites = 50000))

    expect_error(object = imputeLDKNNi(filterGenoObj, knnTaxa = 1))

    expect_error(object = imputeLDKNNi(filterGenoObj, knnTaxa = -1))

    expect_error(object = imputeLDKNNi(filterGenoObj, knnTaxa = 200000))

    expect_error(imputeLDKNNi(filterGenoObj, knnTaxa = "1"))

    expect_error(imputeLDKNNi(filterGenoObj, highLDSSites = "1"))
})


test_that("imputeLDKNNi() returns correct data", {
    expect_true(
        inherits(
            imputeLDKNNi(
                tasObj = filterGenoObj
            ),
            "TasselGenotypePhenotype"
        )
    )
})


