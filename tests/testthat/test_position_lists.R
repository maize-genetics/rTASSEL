# === Tests for position list functions =============================

## Preamble - load data ----

### Start logging info
startLogger()
library(rJava)

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
test_that("getPositionList() throws general exceptions.", {
    tmp <- getPositionList(mtcars)
    expect_true(rJava::is.jnull(tmp))
})



## Return tests ----
test_that("position list methods return correct data and classes.", {
    tmp <- getPositionList(tasGeno@jPositionList)
    expect_true(tmp %instanceof% "net.maizegenetics.dna.map.PositionArrayList")

    tmp <- genomicRanges(tasGeno)
    expect_true(class(tmp)[1] == "GRanges")
    expect_true(length(tmp$tasselIndex) == 3093)

})


