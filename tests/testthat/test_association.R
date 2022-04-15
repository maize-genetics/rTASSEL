# === Tests for `assocModelFitter()` ================================

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

### Create kinship object
tasKin <- rTASSEL::kinshipMatrix(tasGenoPhenoFast)

### Association objects
tasBLUE <- assocModelFitter(
    tasObj  = tasPheno,
    formula = . ~ .
)
tasGLM <- assocModelFitter(
    tasObj     = tasGenoPhenoFast,
    formula    = . ~ .,
    fitMarkers = TRUE
)
tasMLM <- assocModelFitter(
    tasObj     = tasGenoPhenoFast,
    formula    = . ~ .,
    fitMarkers = TRUE,
    kinship    = tasKin
)
tasFast <- assocModelFitter(
    tasObj          = tasGenoPhenoFast,
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
       regexp = "tasObj is not of class \"TasselGenotypePhenotype\""
    )
    expect_error(
        object = assocModelFitter(
            tasObj          = tasGeno,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = TRUE,
            maxP            = 1
        ),
       regexp = "tasObj does not contain a Phenotype object"
    )
    expect_error(
        object = assocModelFitter(
            tasObj          = tasGenoPhenoFast,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = TRUE,
            maxP            = "1"
        ),
        regexp = "p-value must be numeric"
    )
    expect_that( # <- regex expect and actual are the same?
        object = assocModelFitter(
            tasObj          = tasGenoPhenoFast,
            formula         = . ~ .,
            fastAssociation = TRUE,
            fitMarkers      = TRUE,
            maxP            = -2
        ),
        condition = throws_error()
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
            tasObj          = tasGenoPhenoFast,
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
    expect_identical(class(tasBLUE), "list")
    expect_identical(class(tasGLM), "list")
    expect_identical(class(tasMLM), "list")
    expect_identical(class(tasFast), "list")

    expect_equal(dim(tasBLUE$BLUE), c(298, 4))
    expect_equal(dim(tasBLUE$BLUE_ANOVA), c(3, 9))
    expect_equal(dim(tasGLM$GLM_Stats), c(559, 18))
    expect_equal(dim(tasGLM$GLM_Genotypes), c(1337, 7))
    expect_equal(dim(tasMLM$MLM_Effects), c(21615, 7))
    expect_equal(dim(tasMLM$MLM_Stats), c(9282, 18))
    expect_equal(dim(tasFast$FastAssociation), c(640, 7))

    expect_equal(
        object   = colnames(tasBLUE$BLUE),
        expected = c("Taxa", "EarHT", "dpoll", "EarDia")
    )
    expect_equal(
        object   = colnames(tasBLUE$BLUE_ANOVA),
        expected = c(
            "Trait", "F", "p", "taxaDF", "taxaMS", "errorDF",
            "errorMS", "modelDF", "modelMS"
        )
    )
    expect_equal(
        object   = colnames(tasGLM$GLM_Stats),
        expected = c(
            "Trait", "Marker", "Chr", "Pos", "marker_F", "p",
            "marker_Rsq", "add_F", "add_p", "dom_F", "dom_p",
            "marker_df", "marker_MS", "error_df", "error_MS",
            "model_df", "model_MS", "minorObs"
        )
    )
    expect_equal(
        object   = colnames(tasGLM$GLM_Genotypes),
        expected = c(
            "Trait", "Marker", "Chr", "Pos", "Obs", "Allele", "Estimate"
        )
    )
    expect_equal(
        object   = colnames(tasMLM$MLM_Effects),
        expected = c(
            "Trait", "Marker", "Locus", "Site", "Allele", "Effect", "Obs"
        )
    )
    expect_equal(
        object   = colnames(tasMLM$MLM_Stats),
        expected = c(
            "Trait", "Marker", "Chr", "Pos", "df", "F", "p", "add_effect",
            "add_F", "add_p", "dom_effect", "dom_F", "dom_p", "errordf",
            "MarkerR2", "Genetic_Var", "Residual_Var", "-2LnLikelihood"
        )
    )
    expect_equal(
        object   = colnames(tasFast$FastAssociation),
        expected = c(
            "Trait", "Marker", "Chr", "Pos", "df", "r2", "p"
        )
    )

})










