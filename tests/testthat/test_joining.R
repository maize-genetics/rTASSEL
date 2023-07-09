## Intersect and union tests -----
phA <- readPhenotypeFromDataFrame(
    data.frame(
        taxa = c("a", "b", "c", "d"),
        weight = c(120, 150, 100, 70)
    ), "taxa"
)

phB <- readPhenotypeFromDataFrame(
    data.frame(
        taxa = c("a", "b", "c"),
        height = c(12, 15, 10)
    ), "taxa"
)

test_that("Intersect join returns correct values", {

    intersectPheno <- intersectJoin(c(phA, phB))
    expect_equal(getTaxaIDs(intersectPheno), c("a", "b", "c"))

    phenoAttrib <- extractPhenotypeAttDf(intersectPheno@jPhenotypeTable)
    expect_equal(phenoAttrib$traitName, c("Taxa", "weight", "height"))
})

test_that("Union join returns correct values", {
    unionPheno <- unionJoin(c(phA, phB))
    expect_equal(getTaxaIDs(unionPheno), c("a", "b", "c", "d"))

    phenoAttrib <- extractPhenotypeAttDf(unionPheno@jPhenotypeTable)
    expect_equal(phenoAttrib$traitName, c("Taxa", "weight", "height"))
})



## Concatenation tests ----
phA1 <-readPhenotypeFromDataFrame(
    data.frame(
        taxa = c("a", "b", "c"),
        height = c(12, 15, 10)
    ), "taxa"
)

phA2 <- readPhenotypeFromDataFrame(
    data.frame(
        taxa = c("d", "e", "f", "g"),
        height = c(14, 50, 13, 23)
    ), "taxa"
)

test_that("Concatenation returns correct values", {
    concatPheno <- concatenate(c(phA1, phA2))
    expect_equal(getTaxaIDs(concatPheno), c("a", "b", "c", "d", "e", "f", "g"))

    phenoAttrib <- extractPhenotypeAttDf(concatPheno@jPhenotypeTable)
    expect_equal(phenoAttrib$traitName, c("Taxa", "height"))
})



## PCA testing ----
genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)
phenoPath <- system.file(
    "extdata",
    "mdp_traits_nomissing.txt",
    package = "rTASSEL"
)
tasGeno <- readGenotypeTableFromPath(
    path = genoPathHMP
)
tasPheno <- readPhenotypeFromPath(
    path = phenoPath
)
pcaRes <- pca(tasGeno)

test_that("Joining returns correct values with PCA objects", {
    intersectPheno <- intersectJoin(c(pcaRes, tasPheno))
    phenoAttrib <- rTASSEL:::extractPhenotypeAttDf(intersectPheno@jPhenotypeTable)
    expect_equal(
        phenoAttrib$traitName,
        c("Taxa", "EarHT", "dpoll", "EarDia", "PC1", "PC2", "PC3", "PC4", "PC5")
    )

    unionPheno <- unionJoin(c(pcaRes, tasPheno))
    phenoAttrib <- extractPhenotypeAttDf(unionPheno@jPhenotypeTable)
    expect_equal(
        phenoAttrib$traitName,
        c("Taxa", "EarHT", "dpoll", "EarDia", "PC1", "PC2", "PC3", "PC4", "PC5")
    )
})


