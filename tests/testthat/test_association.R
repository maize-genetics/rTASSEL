# === Tests for `assocModelFitter()` ================================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasPheno   <- rtObjs$ph_nomiss
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss

### Create kinship object
tasKin <- rTASSEL::kinshipMatrix(tasDataset)

### Create error-prone kinship object
set.seed(123)
m <- 10
s <- matrix(rnorm(100), m)
s[lower.tri(s)] <- t(s)[lower.tri(s)]
diag(s) <- 2
colnames(s) <- rownames(s) <- paste0("s_", seq_len(m))
tasKinError <- asTasselDistanceMatrix(s)

### Multi data type pheno object (factors + covariates)
tasDatasetCov <- rtObjs$ds_hmp_ph_full

### Create no missing data set with factors/cov (fast association)
fullDf <- as.data.frame(rtObjs$ph_full)
noMissingDf <- fullDf[
    which(!is.na(fullDf$EarHT)),
    c("Taxa", "location", "EarHT", "Q1", "Q2", "Q3")
]
tasDatasetCovNoMiss <- readGenomicDataset(
    genotype  = tasGeno,
    phenotype = noMissingDf,
    attr      = data.frame(
        col_id      = colnames(noMissingDf),
        tassel_attr = c("taxa", "factor", "data", rep("covariate", 3))
    )
)

### Association objects
tasBLUE <- assocModelFitter(
    tasObj  = tasPheno,
    formula = . ~ .
)
tasGLM <- assocModelFitter(
    tasObj     = tasDataset,
    formula    = . ~ .,
    fitMarkers = TRUE
)
tasMLM <- assocModelFitter(
    tasObj     = tasDataset,
    formula    = . ~ .,
    fitMarkers = TRUE,
    kinship    = tasKin
)
tasFast <- assocModelFitter(
    tasObj          = tasDataset,
    formula         = . ~ .,
    fitMarkers      = TRUE,
    fastAssociation = TRUE
)


## Error tests ----

### General errors ----
test_that("assocModelFitter() throws general exceptions.", {
    expect_error(
        object = assocModelFitter(
            tasObj          = mtcars,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = TRUE,
            maxP            = 1
        ),
       regexp = "Unsupported input object"
    )
    expect_error(
        object = assocModelFitter(
            tasObj          = tasGeno,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = TRUE,
            maxP            = 1
        ),
       regexp = "needs phenotype data"
    )
    expect_error(
        object = assocModelFitter(
            tasObj          = tasDataset,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = TRUE,
            maxP            = "1"
        ),
        regexp = "p-value must be numeric"
    )
    expect_that( # <- regex expect and actual are the same?
        object = assocModelFitter(
            tasObj          = tasDataset,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = TRUE,
            maxP            = -2
        ),
        condition = throws_error()
    )
})

test_that("assocModelFitter() kinship parameter throws correct exceptions", {
    expect_error(
        object = assocModelFitter(
            tasObj          = tasDataset,
            formula         = . ~ .,
            kinship         = mtcars,
            fitMarkers      = TRUE,
            maxP            = -2
        )
    )
    expect_error(
        object = assocModelFitter(
            tasObj          = tasDataset,
            formula         = . ~ .,
            kinship         = tasKinError,
            fitMarkers      = TRUE,
            maxP            = -2
        )
    )
})

test_that("assocModelFitter() formula parameter throw correct exceptions", {
    expect_error(
        object = assocModelFitter(
            tasObj     = tasDataset,
            formula    = not_a_trait ~ .,
            fitMarkers = TRUE
        )
    )
    expect_error(
        object = assocModelFitter(
            tasObj     = tasDatasetCov,
            formula    = Q1 ~ .,
            fitMarkers = TRUE
        )
    )
    expect_error(
        object = assocModelFitter(
            tasObj     = tasDatasetCov,
            formula    = . ~ EarHT,
            fitMarkers = TRUE
        )
    )
})


### Specific errors - BLUEs ----
test_that("assocModelFitter() throws exceptions (BLUEs).", {
    expect_that(
        object = assocModelFitter(
            tasObj          = tasGeno,
            formula         = . ~ .,
            fastAssociation = FALSE,
            fitMarkers      = FALSE,
            maxP            = 1
        ),
        condition = throws_error()
    )
    expect_that(
        object = assocModelFitter(
            tasObj          = tasGeno,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = FALSE,
            maxP            = 1
        ),
        condition = throws_error()
    )
    expect_that(
        object = assocModelFitter(
            tasObj          = tasGeno,
            formula         = . ~ .,
            fastAssociation = FALSE,
            fitMarkers      = TRUE,
            maxP            = 1
        ),
        condition = throws_error()
    )
})

### Specific errors - GLM ----
test_that("assocModelFitter() throws exceptions (GLM).", {
    expect_that(
        object = assocModelFitter(
            tasObj          = tasPheno,
            formula         = . ~ .,
            fastAssociation = FALSE,
            fitMarkers      = TRUE,
            maxP            = 1
        ),
        condition = throws_error()
    )
})

### Specific errors - MLM ----
test_that("assocModelFitter() throws exceptions (MLM).", {
    expect_that(
        object = assocModelFitter(
            tasObj          = tasPheno,
            formula         = . ~ .,
            fastAssociation = FALSE,
            fitMarkers      = TRUE,
            kinship         = tasKin,
            maxP            = 1
        ),
        condition = throws_error()
    )
})

### Specific errors - fast association ----
test_that("assocModelFitter() throws exceptions (fast association).", {
    expect_that(
        object = assocModelFitter(
            tasObj          = tasPheno,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = TRUE,
            maxP            = 1
        ),
        condition = throws_error()
    )
    expect_that(
        object = assocModelFitter(
            tasObj          = tasDataset,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = TRUE,
            kinship         = tasKin,
            maxP            = 1
        ),
        condition = throws_error()
    )
})


## Return tests ----

### Return data objects - general ----
test_that("BLUE analysis return correct data types.", {
    expect_true(is(tasBLUE, "AssociationResults"))
    expect_true(is(tasGLM, "AssociationResults"))
    expect_true(is(tasMLM, "AssociationResults"))
    expect_true(is(tasFast, "AssociationResults"))

    # expect_identical(class(tasBLUE), "AssociationResults")
    # expect_identical(class(tasGLM), "AssociationResults")
    # expect_identical(class(tasMLM), "AssociationResults")
    # expect_identical(class(tasFast), "AssociationResults")

    expect_equal(dim(tableReport(tasBLUE, "BLUE")), c(298, 4))
    expect_equal(dim(tableReport(tasBLUE, "BLUE_ANOVA")), c(3, 9))
    expect_equal(dim(tableReport(tasGLM, "GLM_Stats")), c(559, 18))
    expect_equal(dim(tableReport(tasGLM, "GLM_Genotypes")), c(1337, 7))
    expect_equal(dim(tableReport(tasMLM, "MLM_Effects")), c(21615, 7))
    expect_equal(dim(tableReport(tasMLM, "MLM_Stats")), c(9282, 18))
    expect_equal(dim(tableReport(tasFast, "FastAssociation")), c(640, 7))

    expect_equal(
        object   = colnames(tableReport(tasBLUE, "BLUE")),
        expected = c("Taxa", "EarHT", "dpoll", "EarDia")
    )
    expect_equal(
        object   = colnames(tableReport(tasBLUE, "BLUE_ANOVA")),
        expected = c(
            "Trait", "F", "p", "taxaDF", "taxaMS", "errorDF",
            "errorMS", "modelDF", "modelMS"
        )
    )
    expect_equal(
        object   = colnames(tableReport(tasGLM, "GLM_Stats")),
        expected = c(
            "Trait", "Marker", "Chr", "Pos", "marker_F", "p",
            "marker_Rsq", "add_F", "add_p", "dom_F", "dom_p",
            "marker_df", "marker_MS", "error_df", "error_MS",
            "model_df", "model_MS", "minorObs"
        )
    )
    expect_equal(
        object   = colnames(tableReport(tasGLM, "GLM_Genotypes")),
        expected = c(
            "Trait", "Marker", "Chr", "Pos", "Obs", "Allele", "Estimate"
        )
    )
    expect_equal(
        object   = colnames(tableReport(tasMLM, "MLM_Effects")),
        expected = c(
            "Trait", "Marker", "Locus", "Site", "Allele", "Effect", "Obs"
        )
    )
    expect_equal(
        object   = colnames(tableReport(tasMLM, "MLM_Stats")),
        expected = c(
            "Trait", "Marker", "Chr", "Pos", "df", "F", "p", "add_effect",
            "add_F", "add_p", "dom_effect", "dom_F", "dom_p", "errordf",
            "MarkerR2", "Genetic_Var", "Residual_Var", "-2LnLikelihood"
        )
    )
    expect_equal(
        object   = colnames(tableReport(tasFast, "FastAssociation")),
        expected = c(
            "Trait", "Marker", "Chr", "Pos", "df", "r2", "p"
        )
    )
})


## Miscellaneous logic checks ----
test_that("assocModelFitter() handles threads", {
    expect_message(
        assocModelFitter(
            tasObj          = tasDataset,
            formula         = . ~ .,
            fitMarkers      = TRUE,
            fastAssociation = TRUE,
            maxThreads      = 1
        )
    )
})

test_that("assocModelFitter() saves to disk", {
    tmpOut <- paste0(tempdir(), "/test_prefix")
    tmpObj <- assocModelFitter(
        tasObj     = tasDataset,
        formula    = . ~ .,
        fitMarkers = TRUE,
        outputFile = tmpOut
    )
    expect_true(file.exists(paste0(tmpOut, "_allele.txt")))
    expect_true(file.exists(paste0(tmpOut, "_site.txt")))
})

test_that("assocModelFitter() order of operations is correct", {
    expect_message(
        object = assocModelFitter(
            tasObj     = tasDatasetCovNoMiss,
            formula    = . ~ .,
            fitMarkers = TRUE,
            fastAssociation = TRUE
        )
    )
    expect_message(
        object = assocModelFitter(
            tasObj     = tasDatasetCovNoMiss,
            formula    = . ~ .,
            fitMarkers = TRUE,
            kinship = tasKin
        )
    )

    expect_message(
        object = assocModelFitter(
            tasObj     = tasDatasetCovNoMiss,
            formula    = EarHT ~ location,
            fitMarkers = TRUE
        )
    )
})


## Back-compatibility ----
test_that("assocModelFitter() accepts a deprecated TasselGenotypePhenotype", {
    legacyGLM <- assocModelFitter(
        tasObj     = rtObjsLegacy$gt_hmp_ph_nomiss,
        formula    = . ~ .,
        fitMarkers = TRUE
    )

    expect_s4_class(legacyGLM, "AssociationResults")
    expect_equal(
        tableReport(legacyGLM, "GLM_Stats"),
        tableReport(tasGLM, "GLM_Stats")
    )
})


