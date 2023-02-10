# === Tests for visualization methods for "Diversity" analysis ======

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
# TODO - make these not fail on servers...
# test_that("ldJavaApp() throws the right exceptions", {
#     expect_error(
#         object = ldJavaApp(
#             tasObj = mtcars
#         ),
#         regexp = "tasObj is not of class \"TasselGenotypePhenotype\""
#     )
#
#     expect_error(
#         object = ldJavaApp(
#             tasObj = tasPheno
#         ),
#         regexp = "tasObj does contain a Genotype object"
#     )
#
#     expect_message(
#         object = ldJavaApp(
#             tasObj = tasGeno,
#             ldType = "SlidingWindow"
#         ),
#         regexp = "`windowSize` is not set - setting to `1`"
#     )
#
#     expect_message(
#         object = ldJavaApp(
#             tasObj = tasGeno,
#             ldType = "All"
#         ),
#         regexp = "This *might* produce a massive matrix. You have been warned."
#     )
#
# })


test_that("ldPlot() throws the right exceptions", {
    expect_error(
        object = ldPlot(
            tasObj = mtcars
        )
    )
})


test_that("ldPlot() returns correct data types", {
    tasGenoPhenoFilt <- filterGenotypeTableSites(
        tasObj              = tasGeno,
        siteRangeFilterType = "position",
        startPos            = 228e6,
        endPos              = 300e6,
        startChr            = 2,
        endChr              = 2
    )

    myLD <- ldPlot(
        tasObj  = tasGenoPhenoFilt,
        ldType  = "All",
        plotVal = "r2",
        verbose = FALSE
    )
    myLDSW <- ldPlot(
        tasObj  = tasGenoPhenoFilt,
        ldType  = "SlidingWindow",
        plotVal = "r2",
        verbose = FALSE
    )

    truthLabels <- c("x", "y", "val", "group")
    truthObs    <- 544
    truthObsSW  <- 60

    expect_true(inherits(myLD, "ggplot"))
    expect_equal(as.vector(unlist(myLD$labels)), truthLabels)
    expect_equal(nrow(myLD$data), truthObs)
    expect_equal(nrow(myLDSW$data[!is.na(myLDSW$data$val), ]), truthObsSW)
})







