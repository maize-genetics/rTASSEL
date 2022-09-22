
genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)
tasGeno <- readGenotypeTableFromPath(
    path = genoPathHMP
)


test_that("Genotype table coercion to matrix returns correct data", {
    m <- as.matrix(tasGeno)

    expect_equal(dim(m), c(281, 3093))

})


