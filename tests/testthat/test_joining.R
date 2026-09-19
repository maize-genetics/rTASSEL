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


phC <- readPhenotype(
    data.frame(
        taxa = c("a", "b", "c"),
        girth = c(3, 4, 5)
    ),
    attr = phAttr("girth")
)

test_that("Joins take any number of objects without a list", {
    intersectPheno <- intersectJoin(phA, phB, phC)

    expect_s4_class(intersectPheno, "TasselPhenotype")
    expect_equal(getTaxaIDs(intersectPheno), c("a", "b", "c"))
    expect_equal(
        attributeData(intersectPheno)$trait_id,
        c("Taxa", "weight", "height", "girth")
    )

    unionPheno <- unionJoin(phA, phB, phC)
    expect_equal(getTaxaIDs(unionPheno), c("a", "b", "c", "d"))
})

test_that("Variadic and list call styles agree", {
    expect_equal(
        attributeData(intersectJoin(phA, phB)),
        attributeData(intersectJoin(c(phA, phB)))
    )
    expect_equal(
        getTaxaIDs(unionJoin(phA, phB)),
        getTaxaIDs(unionJoin(list(phA, phB)))
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

test_that("Joining returns correct values with MDS objects", {
    mdsRes   <- mds(distanceMatrix(rtObjs$gt_hmp))
    tasPheno <- rtObjs$ph_nomiss

    expectedTraits <- c(
        "Taxa", "PC1", "PC2", "PC3", "PC4", "PC5",
        "EarHT", "dpoll", "EarDia"
    )

    intersectPheno <- intersectJoin(mdsRes, tasPheno)
    expect_s4_class(intersectPheno, "TasselPhenotype")
    expect_equal(attributeData(intersectPheno)$trait_id, expectedTraits)
})

test_that("Joining returns correct values with BLUE objects", {
    blueRes <- assocModelFitter(rtObjs$ph_nomiss, . ~ .)
    phCov   <- readPhenotype(rtFiles$ph_popstruct_path)

    joined <- intersectJoin(blueRes, phCov)

    expect_s4_class(joined, "TasselPhenotype")
    expect_equal(
        attributeData(joined)$trait_id,
        c("Taxa", "EarHT", "dpoll", "EarDia", "Q1", "Q2", "Q3")
    )
    expect_equal(
        attributeData(joined)$trait_type,
        c("taxa", rep("data", 3), rep("covariate", 3))
    )

    # The estimates themselves, rather than a table of the same shape
    blueDf   <- tableReport(blueRes, "BLUE")
    joinedDf <- as.data.frame(joined)
    expect_equal(
        joinedDf$EarHT,
        blueDf$EarHT[match(joinedDf$Taxa, blueDf$Taxa)]
    )

    withGt <- intersectJoin(rtObjs$gt_hmp, blueRes, phCov)
    expect_s4_class(withGt, "TasselGenomicDataset")
    expect_equal(
        traitNames(withGt),
        c("EarHT", "dpoll", "EarDia", "Q1", "Q2", "Q3")
    )
})

test_that("Joins reject association results without phenotype data", {
    mockResults <- list("td_1" = iris)

    expect_error(
        intersectJoin(
            methods::new(
                "AssociationResultsGLM",
                results   = mockResults,
                traits    = "trait_1",
                assocType = "GLM"
            ),
            phA
        ),
        "cannot join <AssociationResultsGLM> results"
    )

    # A BLUE object built by hand has no TASSEL phenotype behind it
    expect_error(
        intersectJoin(
            methods::new(
                "AssociationResultsBLUE",
                results   = mockResults,
                traits    = "trait_1",
                assocType = "BLUE"
            ),
            phA
        ),
        "no phenotype data"
    )
})

test_that("Joining accepts a genomic dataset's phenotype data", {
    joined <- intersectJoin(c(rtObjs$ds_hmp_ph_nomiss, pca(rtObjs$gt_hmp)))

    expect_s4_class(joined, "TasselPhenotype")
    expect_equal(
        attributeData(joined)$trait_id,
        c("Taxa", "EarHT", "dpoll", "EarDia", paste0("PC", 1:5))
    )
})

test_that("Joining a genotype object returns a genomic dataset", {
    phCov <- readPhenotype(rtFiles$ph_popstruct_path)

    joined <- intersectJoin(rtObjs$gt_hmp, rtObjs$ph_nomiss, phCov)

    expect_s4_class(joined, "TasselGenomicDataset")
    expect_equal(
        traitNames(joined),
        c("EarHT", "dpoll", "EarDia", "Q1", "Q2", "Q3")
    )
    expect_equal(
        joined@genotype@jRefObj$numberOfSites(),
        rtObjs$gt_hmp@jRefObj$numberOfSites()
    )
    expect_setequal(
        getTaxaIDs(joined),
        intersect(getTaxaIDs(rtObjs$ds_hmp_ph_nomiss), getTaxaIDs(phCov))
    )
})

test_that("The join mode carries through to the genotype join", {
    # A taxon the genotype table does not know about, so the genotype
    # join mode - not just the phenotype one - decides whether it stays
    phTaxa <- c(head(getTaxaIDs(rtObjs$gt_hmp), 3), "fake_line")

    phX <- readPhenotype(
        data.frame(taxa = phTaxa, weight = c(120, 150, 100, 70)),
        attr = phAttr("weight")
    )
    phY <- readPhenotype(
        data.frame(taxa = phTaxa, height = c(12, 15, 10, 9)),
        attr = phAttr("height")
    )

    intersectDs <- intersectJoin(rtObjs$gt_hmp, phX, phY)
    unionDs     <- unionJoin(rtObjs$gt_hmp, phX, phY)

    expect_s4_class(unionDs, "TasselGenomicDataset")
    expect_equal(getTaxaIDs(phenotype(intersectDs)), head(phTaxa, 3))
    expect_equal(getTaxaIDs(phenotype(unionDs)), phTaxa)
})

test_that("Joins reject empty and unsupported input", {
    expect_error(intersectJoin(list()), "at least one object")
    expect_error(intersectJoin(c(phA, mtcars)), "Unsupported input object")
})

test_that("Joins reject genotype input they cannot use", {
    expect_error(
        concatenate(rtObjs$gt_hmp, phA1, phA2),
        "does not accept genotype data"
    )
    expect_error(
        intersectJoin(rtObjs$gt_hmp, readGenotype(rtMatrices$num_gt_sm)),
        "Only genotype tables holding allele calls"
    )
})


## Genotype table joining ----
gtChr1 <- readGenotype(returnSysFiles("rt_sub_chr1.vcf"))
gtChr5 <- readGenotype(returnSysFiles("rt_sub_chr5.vcf"))

test_that("Joining genotype tables returns their sites in genomic order", {
    joined <- intersectJoin(gtChr1, gtChr5)

    expect_s4_class(joined, "TasselGenotype")
    expect_equal(getTaxaIDs(joined), getTaxaIDs(gtChr1))
    expect_equal(
        as.data.frame(positionList(joined))[, c("Name", "Chromosome")],
        rbind(
            as.data.frame(positionList(gtChr1)),
            as.data.frame(positionList(gtChr5))
        )[, c("Name", "Chromosome")]
    )
    expect_equal(
        as.matrix(joined),
        cbind(as.matrix(gtChr1), as.matrix(gtChr5))
    )

    # The sites TASSEL reports are only as good as the site scores it can
    # compute from them, which a lazy view over the tables cannot
    expect_equal(dim(as.matrix(kinshipMatrix(joined))), c(5L, 5L))
})

test_that("The order genotype tables are given in does not matter", {
    expect_equal(
        as.matrix(intersectJoin(gtChr5, gtChr1)),
        as.matrix(intersectJoin(gtChr1, gtChr5))
    )
})

test_that("Joining genotype tables handles interleaved sites", {
    nSites <- nrow(positionList(rtObjs$gt_hmp))

    oddSites  <- rtObjs$gt_hmp[, sites(seq(1, nSites, by = 2))]
    evenSites <- rtObjs$gt_hmp[, sites(seq(2, nSites, by = 2))]

    rejoined <- intersectJoin(oddSites, evenSites)

    expect_equal(
        positionList(rejoined),
        positionList(rtObjs$gt_hmp)
    )
    expect_equal(as.matrix(rejoined), as.matrix(rtObjs$gt_hmp))
})

test_that("The join mode decides which taxa a genotype join keeps", {
    gtChr5Sub <- gtChr5[taxa("33-16", "38-11"), ]

    intersectGt <- intersectJoin(gtChr1, gtChr5Sub)
    unionGt     <- unionJoin(gtChr1, gtChr5Sub)

    expect_equal(getTaxaIDs(intersectGt), c("33-16", "38-11"))
    expect_equal(getTaxaIDs(unionGt), getTaxaIDs(gtChr1))

    # The calls the subset never made come back as missing
    unionCalls <- as.matrix(unionGt, type = "allele")
    chr5Sites  <- nrow(positionList(gtChr1)) + seq_len(nrow(positionList(gtChr5)))
    expect_true(all(unionCalls[c("4226", "4722", "A188"), chr5Sites] == "N"))
    expect_equal(
        unionCalls[c("33-16", "38-11"), chr5Sites],
        as.matrix(gtChr5Sub, type = "allele")
    )
})

test_that("Genotype tables and phenotype data can be joined at once", {
    joined <- intersectJoin(gtChr1, gtChr5, rtObjs$ph_nomiss)

    expect_s4_class(joined, "TasselGenomicDataset")
    expect_equal(nrow(positionList(joined)), 17)
    expect_equal(traitNames(joined), c("EarHT", "dpoll", "EarDia"))
    expect_setequal(
        getTaxaIDs(joined),
        intersect(getTaxaIDs(gtChr1), getTaxaIDs(rtObjs$ph_nomiss))
    )
})

test_that("A single genotype object is returned as it was given", {
    expect_equal(as.matrix(intersectJoin(gtChr1)), as.matrix(gtChr1))
    expect_s4_class(intersectJoin(gtChr1), "TasselGenotype")
})

test_that("Genotype tables sharing no taxa cannot be joined", {
    expect_error(
        intersectJoin(gtChr1[taxa("33-16"), ], gtChr5[taxa("4226"), ]),
        "Could not join the genotype tables"
    )
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

    legacyGtJoin <- intersectJoin(
        rtObjsLegacy$gt_hmp,
        rtObjsLegacy$ph_nomiss,
        readPhenotypeFromPath(rtFiles$ph_popstruct_path)
    )
    expect_s4_class(legacyGtJoin, "TasselGenotypePhenotype")
    expect_false(rJava::is.jnull(getGenotypeTable(legacyGtJoin)))

    legacyGtOnlyJoin <- intersectJoin(
        readGenotypeTableFromPath(returnSysFiles("rt_sub_chr1.vcf")),
        readGenotypeTableFromPath(returnSysFiles("rt_sub_chr5.vcf"))
    )
    expect_s4_class(legacyGtOnlyJoin, "TasselGenotypePhenotype")
    expect_equal(nrow(positionList(legacyGtOnlyJoin)), 17)

    legacyMerge <- mergeGenotypeTables(list(
        rtObjsLegacy$gt_hmp,
        rtObjsLegacy$gt_hmp
    ))
    expect_s4_class(legacyMerge, "TasselGenotypePhenotype")
})
