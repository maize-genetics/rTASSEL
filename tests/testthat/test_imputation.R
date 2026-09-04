# === Tests for imputation methods ==================================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasPheno   <- rtObjs$ph_nomiss
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss

### Imputation is slow, so cut both fixtures down to a small block of
### sites and the taxa whose IDs start with a digit or an "A"
smallTaxa  <- taxaWhere(grepl("^[0-9]|^A", taxaId))
smallSites <- sites(1:11)

filterGenoObj      <- tasGeno[smallTaxa, smallSites]
filterGenoPhenoObj <- tasDataset[smallTaxa, smallSites]


## Imputation (Numeric) ----
test_that("imputeNumeric() returns correct exceptions", {
    expect_error(
        object = imputeNumeric(mtcars),
        regexp = "Unsupported input object"
    )

    expect_error(
        object = imputeNumeric(tasPheno),
        regexp = "needs genotype data"
    )
})


test_that("imputeNumeric() returns a genotype for every distance metric", {
    expect_s4_class(
        imputeNumeric(
            tasObj = filterGenoObj,
            byMean = TRUE,
            nearestNeighbors = 5,
            distance = "Euclidean"
        ),
        "TasselGenotype"
    )

    for (dist in c("Euclidean", "Manhattan", "Cosine")) {
        expect_s4_class(
            imputeNumeric(
                tasObj = filterGenoObj,
                byMean = FALSE,
                nearestNeighbors = 5,
                distance = dist
            ),
            "TasselGenotype"
        )
    }
})


test_that("imputeNumeric() returns a dataset for every distance metric", {
    expect_s4_class(
        imputeNumeric(
            tasObj = filterGenoPhenoObj,
            byMean = TRUE,
            nearestNeighbors = 5,
            distance = "Euclidean"
        ),
        "TasselGenomicDataset"
    )

    for (dist in c("Euclidean", "Manhattan", "Cosine")) {
        expect_s4_class(
            imputeNumeric(
                tasObj = filterGenoPhenoObj,
                byMean = FALSE,
                nearestNeighbors = 5,
                distance = dist
            ),
            "TasselGenomicDataset"
        )
    }
})


test_that("imputeNumeric() keeps a dataset's phenotype data attached", {
    imputed <- imputeNumeric(filterGenoPhenoObj)

    expect_equal(traitNames(imputed), traitNames(filterGenoPhenoObj))
    expect_equal(
        nrow(as.data.frame(imputed)),
        nrow(as.data.frame(filterGenoPhenoObj))
    )
})



## Imputation (LD KNNi) ----
test_that("imputeLDKNNi() returns correct exceptions", {
    expect_error(
        object = imputeLDKNNi(mtcars),
        regexp = "Unsupported input object"
    )

    expect_error(
        object = imputeLDKNNi(tasPheno),
        regexp = "needs genotype data"
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
    expect_s4_class(imputeLDKNNi(filterGenoObj), "TasselGenotype")
    expect_s4_class(imputeLDKNNi(filterGenoPhenoObj), "TasselGenomicDataset")
})


## Back-compatibility ----
test_that("imputation accepts a deprecated TasselGenotypePhenotype", {
    legacySub <- filterGenotypeTableSites(
        tasObj = rtObjsLegacy$gt_hmp_ph_nomiss,
        siteRangeFilterType = "sites",
        startSite = 0,
        endSite = 10
    )

    expect_s4_class(imputeNumeric(legacySub), "TasselGenotypePhenotype")
    expect_s4_class(imputeLDKNNi(legacySub), "TasselGenotypePhenotype")
})


