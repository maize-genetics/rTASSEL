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
})









