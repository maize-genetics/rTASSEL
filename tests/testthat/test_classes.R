# === Tests for rTASSEL classes =====================================

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



## Return tests ----
test_that("rTASSEL read functions return correct data type.", {
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

    expect_s4_class(
        object = tasPheno,
        class  = "TasselGenotypePhenotype"
    )
    expect_s4_class(
        object = tasGeno,
        class  = "TasselGenotypePhenotype"
    )
    expect_s4_class(
        object = tasGenoPheno,
        class  = "TasselGenotypePhenotype"
    )

    tmp <- methods::getSlots("TasselGenotypePhenotype")
    expect_equal(
        object = names(tmp),
        expected = c(
            "name", "jTasselObj", "jTaxaList", "jPositionList",
            "jGenotypeTable", "jPhenotypeTable"
        )
    )
    tmp <- as.vector(tmp)
    expect_equal(
        object = tmp,
        expected = c("character", rep("jobjRef", 5))
    )

    showGeno <- utils::capture.output(tasGeno)
    showPheno <- utils::capture.output(tasPheno)
    expect_equal(length(showGeno), 8)
    expect_equal(length(showPheno), 10)
    expect_equal(showGeno[7], "  Genotype Table..... [x]")
    expect_equal(showGeno[8], "  Phenotype Table.... [ ]")
    expect_equal(showPheno[7], "  Genotype Table..... [ ]")
    expect_equal(showPheno[8], "  Phenotype Table.... [x]")

})

test_that("TasselDistanceMatrix object return correct output", {
    tasGeno <- readGenotypeTableFromPath(
        path = genoPathHMP
    )
    tasKin <- kinshipMatrix(tasGeno)

    showKin <- utils::capture.output(tasKin)
    expect_equal(length(showKin), 9)
    expect_equal(showKin[1], "A TasselDistanceMatrix Object of 281 x 281 elements:")
    expect_equal(showKin[2], "")
    expect_equal(showKin[3], "                  33-16      38-11       4226       4722    ...    YU796NS")
    expect_equal(showKin[4], "       33-16     1.7978     0.0396     0.0614    -0.0083    ...    -0.0003")
    expect_equal(showKin[5], "       38-11     0.0396     1.9100     0.0198    -0.0063    ...    -0.0190")
    expect_equal(showKin[6], "        4226     0.0614     0.0198     1.9277    -0.0165    ...     0.1306")
    expect_equal(showKin[7], "        4722    -0.0083    -0.0063    -0.0165     1.4539    ...     0.0331")
    expect_equal(showKin[8], "         ...        ...        ...        ...        ...    ...        ...")
    expect_equal(showKin[9], "     YU796NS    -0.0003    -0.0190     0.1306     0.0331    ...     1.8536")
    expect_equal(length(rownames(tasKin)), 281)
})


test_that("TasselGenotypePhenotype() returns correct NAs", {
    testDf <- data.frame(
        Taxa = c("t_a", "t_b", "t_c"),
        trait = c(0.5, 0.1, 0.9)
    )

    testData <- readPhenotypeFromDataFrame(testDf, "Taxa")
    testData@jTaxaList <- rJava::.jnull()

    testOut <- capture.output(testData)

    expect_true(grepl("NA", testOut[3]))

    set.seed(123)
    testDf <- data.frame(
        Taxa = c("t_a", "t_b", "t_c"),
        trait_1 = rnorm(3),
        trait_2 = rnorm(3),
        trait_3 = rnorm(3),
        trait_4 = rnorm(3),
        trait_5 = rnorm(3),
        trait_6 = rnorm(3),
        trait_7 = rnorm(3)
    )
    testData <- readPhenotypeFromDataFrame(testDf, "Taxa")
    testOut <- capture.output(testData)
    expect_true(grepl("with 3 more IDs", testOut[10]))
})


test_that(".getTASSELClass() returns correct data and exceptions", {
    expect_error(.getTASSELClass(rJava::.jnull(), "non_entity"))
    expect_error(.getTASSELClass(rJava::.jnull(), "GenotypePhenotype"))
    expect_error(.getTASSELClass(rJava::.jnull(), "GenotypeTable"))
    expect_error(.getTASSELClass(rJava::.jnull(), "Phenotype"))
    expect_error(.getTASSELClass(rJava::.jnull(), "TaxaList"))
    expect_error(.getTASSELClass(rJava::.jnull(), "PositionList"))


})
























