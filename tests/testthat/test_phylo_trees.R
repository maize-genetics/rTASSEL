# === Tests for phylogenetic analyses ===============================

## Preamble - load data ----

### Start logging info
startLogger()


### Create rTASSEL genotype only object
tasGeno <- readGenotypeTableFromPath(
    path = system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
)

tasGenoSub <- filterGenotypeTableTaxa(
    tasObj = tasGeno,
    taxa = c("33-16", "38-11", "4226", "4722", "A188", "A214N")
)


## Error tests ----
test_that("createTree() throws general exceptions.", {
    expect_error(
        object = createTree(tasObj = mtcars),
        regexp = "tasObj is not of class \"TasselGenotypePhenotype\""
    )
    expect_error(
        object = createTree(tasObj = tasGeno, clustMethod = "NA"),
        regexp = NULL
    )
})


## Return tests ----
test_that("createTree() returns correct data types", {
    t <- createTree(tasObj = tasGeno)
    tSub <- createTree(tasObj = tasGenoSub)

    expect_s3_class(
        object = t,
        class = "phylo"
    )

    expect_equal(
        object = t$Nnode,
        expected = 279
    )

    expect_equal(
        object = tSub$Nnode,
        expected = 4
    )

    expect_equal(
        object = length(t$tip.label),
        expected = 281
    )

    expect_equal(
        object = length(tSub$tip.label),
        expected = 6
    )
})


