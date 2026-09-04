# === Tests for join methods ========================================

## Preamble - load data ----
startLogger()

# Small phenotypes built straight from data frames, so the joins have
# distinct traits to work with
phAttr <- function(trait) {
    data.frame(
        col_id      = c("taxa", trait),
        tassel_attr = c("taxa", "data")
    )
}

phA <- readPhenotype(
    data.frame(
        taxa = c("a", "b", "c", "d"),
        weight = c(120, 150, 100, 70)
    ),
    attr = phAttr("weight")
)

phB <- readPhenotype(
    data.frame(
        taxa = c("a", "b", "c"),
        height = c(12, 15, 10)
    ),
    attr = phAttr("height")
)


## Intersect and union tests -----
test_that("Intersect join returns correct values", {
    intersectPheno <- intersectJoin(c(phA, phB))

    expect_s4_class(intersectPheno, "TasselPhenotype")
    expect_equal(getTaxaIDs(intersectPheno), c("a", "b", "c"))
    expect_equal(
        attributeData(intersectPheno)$trait_id,
        c("Taxa", "weight", "height")
    )
})

test_that("Union join returns correct values", {
    unionPheno <- unionJoin(c(phA, phB))

    expect_s4_class(unionPheno, "TasselPhenotype")
    expect_equal(getTaxaIDs(unionPheno), c("a", "b", "c", "d"))
    expect_equal(
        attributeData(unionPheno)$trait_id,
        c("Taxa", "weight", "height")
    )
})


## Concatenation tests ----
phA1 <- readPhenotype(
    data.frame(
        taxa = c("a", "b", "c"),
        height = c(12, 15, 10)
    ),
    attr = phAttr("height")
)

phA2 <- readPhenotype(
    data.frame(
        taxa = c("d", "e", "f", "g"),
        height = c(14, 50, 13, 23)
    ),
    attr = phAttr("height")
)

test_that("Concatenation returns correct values", {
    concatPheno <- concatenate(c(phA1, phA2))

    expect_s4_class(concatPheno, "TasselPhenotype")
    expect_equal(getTaxaIDs(concatPheno), c("a", "b", "c", "d", "e", "f", "g"))
    expect_equal(attributeData(concatPheno)$trait_id, c("Taxa", "height"))
})


## Joins with other rTASSEL classes ----
test_that("Joining returns correct values with PCA objects", {
    pcaRes   <- pca(rtObjs$gt_hmp)
    tasPheno <- rtObjs$ph_nomiss

    expectedTraits <- c(
        "Taxa", "PC1", "PC2", "PC3", "PC4", "PC5",
        "EarHT", "dpoll", "EarDia"
    )

    intersectPheno <- intersectJoin(c(pcaRes, tasPheno))
    expect_equal(attributeData(intersectPheno)$trait_id, expectedTraits)

    unionPheno <- unionJoin(c(pcaRes, tasPheno))
    expect_equal(attributeData(unionPheno)$trait_id, expectedTraits)
})

test_that("Joining accepts a genomic dataset's phenotype data", {
    joined <- intersectJoin(c(rtObjs$ds_hmp_ph_nomiss, pca(rtObjs$gt_hmp)))

    expect_s4_class(joined, "TasselPhenotype")
    expect_equal(
        attributeData(joined)$trait_id,
        c("Taxa", "EarHT", "dpoll", "EarDia", paste0("PC", 1:5))
    )
})

test_that("Joins reject empty and unsupported input", {
    expect_error(intersectJoin(list()), "at least one object")
    expect_error(intersectJoin(c(phA, mtcars)), "Unsupported input object")
})


## Genotype table merging ----
test_that("mergeGenotypeTables() tests", {
    gtA <- readGenotype(returnSysFiles("rt_sub_chr1.vcf"))
    gtB <- readGenotype(returnSysFiles("rt_sub_chr5.vcf"))

    gtBFilter <- gtB[taxa("33-16", "38-11"), ]

    gtMerged <- mergeGenotypeTables(list(gtA, gtB))
    gtMergedFilter <- mergeGenotypeTables(list(gtA, gtBFilter))

    expect_s4_class(gtMerged, "TasselGenotype")
    expect_error(mergeGenotypeTables(list(gtA, mtcars)))
    expect_error(mergeGenotypeTables(LETTERS))
    expect_error(mergeGenotypeTables(list()), "at least one object")
    expect_equal(length(taxaList(gtMerged)), 5)
    expect_equal(length(taxaList(gtMergedFilter)), 5)
    expect_equal(nrow(positionList(gtMerged)), 17)
})


## Back-compatibility ----
test_that("joins and merges accept deprecated TasselGenotypePhenotype input", {
    legacyPhA <- readPhenotypeFromDataFrame(
        data.frame(taxa = c("a", "b", "c", "d"), weight = c(120, 150, 100, 70)),
        "taxa"
    )
    legacyPhB <- readPhenotypeFromDataFrame(
        data.frame(taxa = c("a", "b", "c"), height = c(12, 15, 10)),
        "taxa"
    )

    legacyJoin <- intersectJoin(c(legacyPhA, legacyPhB))
    expect_s4_class(legacyJoin, "TasselGenotypePhenotype")
    expect_equal(getTaxaIDs(legacyJoin), c("a", "b", "c"))

    legacyMerge <- mergeGenotypeTables(list(
        rtObjsLegacy$gt_hmp,
        rtObjsLegacy$gt_hmp
    ))
    expect_s4_class(legacyMerge, "TasselGenotypePhenotype")
})
