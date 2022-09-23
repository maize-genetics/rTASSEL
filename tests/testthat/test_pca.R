
genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)
tasGeno <- readGenotypeTableFromPath(
    path = genoPathHMP
)


test_that("pca() returns correct data", {
    pcaRes <- pca(tasGeno)

    expect_equal(
        names(pcaRes),
        c("PC_Datum", "Eigenvalues_Datum", "Eigenvectors_Datum", "jPcaObj")
    )

    expect_equal(
        rJava::.jclass(pcaRes$jPcaObj),
        "net.maizegenetics.phenotype.CorePhenotype"
    )
})






