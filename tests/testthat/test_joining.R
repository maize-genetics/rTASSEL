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
    phenoAttrib <- extractPhenotypeAttDf(intersectPheno@jPhenotypeTable)
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


test_that("mergeGenotypeTables() tests", {
    gtA <- readGenotypeTableFromPath(system.file(
        "extdata",
        "rt_sub_chr1.vcf",
        package = "rTASSEL"
    ))
    gtB <- readGenotypeTableFromPath(system.file(
        "extdata",
        "rt_sub_chr5.vcf",
        package = "rTASSEL"
    ))

    gtBFilter <- filterGenotypeTableTaxa(gtB, taxa = c("33-16", "38-11"))

    gtMerged <- mergeGenotypeTables(list(gtA, gtB))
    gtMergedFilter <- mergeGenotypeTables(list(gtA, gtBFilter))

    expect_true(is(gtMerged, "TasselGenotypePhenotype"))
    expect_error(mergeGenotypeTables(list(gtA, mtcars)))
    expect_error(mergeGenotypeTables(LETTERS))
    expect_equal(gtMerged@jTaxaList$numberOfTaxa(), 5)
    expect_equal(gtMergedFilter@jTaxaList$numberOfTaxa(), 5)
    expect_equal(gtMerged@jPositionList$numberOfSites(), 17)
})


# === Tests for taxa mismatch error handling ========================

test_that("intersectJoin() provides informative error when no taxa match", {
    # Create phenotypes with completely different taxa
    phNoMatch1 <- readPhenotypeFromDataFrame(
        data.frame(
            taxa = c("apple", "banana", "cherry"),
            weight = c(100, 200, 150)
        ), "taxa"
    )

    phNoMatch2 <- readPhenotypeFromDataFrame(
        data.frame(
            taxa = c("dog", "cat", "bird"),
            height = c(50, 30, 10)
        ), "taxa"
    )

    expect_error(
        intersectJoin(c(phNoMatch1, phNoMatch2)),
        regexp = "No common taxa were found"
    )

    # Verify error message contains sample taxa names
    expect_error(
        intersectJoin(c(phNoMatch1, phNoMatch2)),
        regexp = "apple|banana|cherry"
    )
})

test_that("readGenotypePhenotype() provides informative error when taxa don't match", {
    # Create a phenotype with taxa that don't match the VCF file
    # The VCF file has taxa like "M0297:C05F2ACXX:5:250021042"
    phenoNoMatch <- data.frame(
        Taxa = c("nonexistent_taxa_1", "nonexistent_taxa_2", "nonexistent_taxa_3"),
        trait1 = c(1.5, 2.5, 3.5)
    )

    genoPathVCF <- system.file(
        "extdata",
        "maize_chr9_10thin40000.recode.vcf",
        package = "rTASSEL"
    )

    expect_error(
        readGenotypePhenotype(
            genoPathOrObj = genoPathVCF,
            phenoPathDFOrObj = phenoNoMatch,
            taxaID = "Taxa"
        ),
        regexp = "No common taxa found"
    )

    # Verify error contains helpful information about sample taxa
    expect_error(
        readGenotypePhenotype(
            genoPathOrObj = genoPathVCF,
            phenoPathDFOrObj = phenoNoMatch,
            taxaID = "Taxa"
        ),
        regexp = "Genotype taxa|Phenotype taxa"
    )
})

test_that("combineTasselGenotypePhenotype() error message shows sample taxa", {
    # Create genotype and phenotype objects with non-matching taxa
    genoPath <- system.file(
        "extdata",
        "maize_chr9_10thin40000.recode.vcf",
        package = "rTASSEL"
    )
    tasGenoVCF <- readGenotypeTableFromPath(genoPath)

    phenoNoMatch <- readPhenotypeFromDataFrame(
        data.frame(
            taxa = c("fake_taxa_1", "fake_taxa_2"),
            value = c(10, 20)
        ), "taxa"
    )

    # Get the Java objects
    genoTable <- getGenotypeTable(tasGenoVCF)
    phenoTable <- getPhenotypeTable(phenoNoMatch)

    # Test that the error includes taxa samples
    expect_error(
        combineTasselGenotypePhenotype(genoTable, phenoTable),
        regexp = "Sample taxa names from each dataset"
    )

    # Verify it mentions case sensitivity
    expect_error(
        combineTasselGenotypePhenotype(genoTable, phenoTable),
        regexp = "case-sensitive"
    )
})

test_that("combineTasselGenotypePhenotype() catches Java IndexOutOfBoundsException", {
    # This test verifies that the tryCatch mechanism properly catches

    # the Java exception and converts it to a helpful R error
    genoPath <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(genoPath)

    # Create phenotype with completely non-matching taxa
    phenoNoMatch <- readPhenotypeFromDataFrame(
        data.frame(
            taxa = c("ZZZZZ_1", "ZZZZZ_2", "ZZZZZ_3"),
            trait = c(1, 2, 3)
        ), "taxa"
    )

    genoTable <- getGenotypeTable(tasGenoHMP)
    phenoTable <- getPhenotypeTable(phenoNoMatch)

    # The error should NOT be the raw Java exception
    expect_error(
        combineTasselGenotypePhenotype(genoTable, phenoTable),
        regexp = "No common taxa found"
    )

    # Verify it does NOT show the raw Java error message
    err <- tryCatch(
        combineTasselGenotypePhenotype(genoTable, phenoTable),
        error = function(e) e
    )
    expect_false(grepl("java.lang.IndexOutOfBoundsException", err$message))
})

test_that("combineTasselGenotypePhenotype() still works when taxa match", {
    # Verify that valid inputs still work correctly after our changes
    genoPath <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    phenoPath <- system.file(
        "extdata",
        "mdp_traits_nomissing.txt",
        package = "rTASSEL"
    )

    tasGenoHMP <- readGenotypeTableFromPath(genoPath)
    tasPhenoValid <- readPhenotypeFromPath(phenoPath)

    genoTable <- getGenotypeTable(tasGenoHMP)
    phenoTable <- getPhenotypeTable(tasPhenoValid)

    # Should succeed without error
    result <- combineTasselGenotypePhenotype(genoTable, phenoTable)

    expect_false(is.null(result))
    expect_false(rJava::is.jnull(result))
})

test_that("readGenotypePhenotype() works with matching taxa from file paths", {
    # Test that the function still works correctly with valid matching data
    genoPath <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    phenoPath <- system.file(
        "extdata",
        "mdp_traits_nomissing.txt",
        package = "rTASSEL"
    )

    # Should succeed without error
    result <- readGenotypePhenotype(
        genoPathOrObj = genoPath,
        phenoPathDFOrObj = phenoPath
    )

    expect_s4_class(result, "TasselGenotypePhenotype")
    expect_true(result@jTaxaList$numberOfTaxa() > 0)
})

test_that("readGenotypePhenotype() error includes actual taxa from VCF file", {
    # Verify that the error message shows the actual taxa names from the VCF
    genoPathVCF <- system.file(
        "extdata",
        "maize_chr9_10thin40000.recode.vcf",
        package = "rTASSEL"
    )

    phenoNoMatch <- data.frame(
        Taxa = c("fake_1", "fake_2"),
        value = c(1, 2)
    )

    err <- tryCatch(
        readGenotypePhenotype(
            genoPathOrObj = genoPathVCF,
            phenoPathDFOrObj = phenoNoMatch,
            taxaID = "Taxa"
        ),
        error = function(e) e
    )

    # Error should contain sample taxa from the VCF file
    # VCF has taxa like "M0297:C05F2ACXX:5:250021042"
    expect_true(grepl("C05F2ACXX|C08L7ACXX", err$message))

    # Error should contain the fake phenotype taxa
    expect_true(grepl("fake_1|fake_2", err$message))
})

test_that("Error handling preserves non-taxa-mismatch Java errors", {
    # Verify that other Java errors are still propagated correctly
    # and not incorrectly caught by our tryCatch

    # Test with NULL genotype table - should give a different error
    phenoValid <- readPhenotypeFromDataFrame(
        data.frame(
            taxa = c("a", "b"),
            value = c(1, 2)
        ), "taxa"
    )

    expect_error(
        combineTasselGenotypePhenotype(NULL, getPhenotypeTable(phenoValid))
    )
})







