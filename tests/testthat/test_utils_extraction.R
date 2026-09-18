# === Tests for JVM -> R extraction helpers and accessors ===========

## Preamble - load data ----

tasGeno    <- rtObjs$gt_hmp
tasPheno   <- rtObjs$ph_full
tasDataset <- rtObjs$ds_hmp_ph_nomiss


## Dosage matrix ----

test_that(".dosageMatrix() labels both axes and maps the missing code", {
    m <- .dosageMatrix(
        javaRefObj(tasGeno),
        taxa      = taxaList(tasGeno),
        siteNames = positionList(tasGeno)$Name
    )

    expect_equal(dim(m), c(281, 3093))
    expect_type(m, "integer")
    expect_equal(rownames(m), taxaList(tasGeno))
    expect_equal(colnames(m), positionList(tasGeno)$Name)

    # TASSEL's unsigned missing value must not survive as a dosage
    expect_false(any(m == 128, na.rm = TRUE))
    expect_true(any(is.na(m)))
    expect_setequal(setdiff(unique(as.vector(m)), NA), c(0L, 1L, 2L))
})

test_that(".dosageMatrix() leaves the raw bytes alone when asked", {
    m <- .dosageMatrix(javaRefObj(tasGeno), asInteger = FALSE)

    expect_type(m, "raw")
    expect_equal(dim(m), c(281, 3093))
    expect_null(dimnames(m))
})

test_that("every dosage caller agrees", {
    expect_equal(as.matrix(tasGeno), as.matrix(rtObjsLegacy$gt_hmp))

    # 'getSumExpFromGenotypeTable()' transposes, since a
    # SummarizedExperiment puts features in rows
    se <- getSumExpFromGenotypeTable(tasGeno, verbose = FALSE)
    expect_equal(
        unname(SummarizedExperiment::assay(se)),
        unname(t(as.matrix(tasGeno)))
    )
})


## Allele string matrix ----

test_that("as.matrix(type = 'allele') matches TASSEL cell for cell", {
    sub <- tasGeno[, sites(1:40)]
    am  <- as.matrix(sub, type = "allele")
    jGt <- javaRefObj(sub)

    expect_type(am, "character")
    expect_equal(dim(am), c(281, 40))
    expect_equal(rownames(am), taxaList(sub))
    expect_equal(colnames(am), positionList(sub)$Name)

    set.seed(42)
    taxonIdx <- sample(nrow(am), 150, replace = TRUE)
    siteIdx  <- sample(ncol(am), 150, replace = TRUE)
    truth <- vapply(
        seq_along(taxonIdx),
        function(k) {
            jGt$genotypeAsString(
                as.integer(taxonIdx[k] - 1L),
                as.integer(siteIdx[k] - 1L)
            )
        },
        FUN.VALUE = character(1)
    )

    expect_equal(am[cbind(taxonIdx, siteIdx)], truth)
})

test_that("the two readings of a call agree on where the data is missing", {
    sub <- tasGeno[, sites(1:200)]

    expect_equal(
        is.na(as.matrix(sub)),
        as.matrix(sub, type = "allele") == "N"
    )
})

test_that("allele coercion holds up on a single taxon and a single site", {
    oneTaxon <- tasGeno[taxaList(tasGeno)[1], sites(1:10)]
    oneSite  <- tasGeno[, sites(5L)]

    expect_equal(dim(as.matrix(oneTaxon, type = "allele")), c(1, 10))
    expect_equal(dim(as.matrix(oneSite, type = "allele")), c(281, 1))

    expect_equal(
        as.vector(as.matrix(oneSite, type = "allele")),
        vapply(
            seq_len(281),
            function(i) javaRefObj(oneSite)$genotypeAsString(as.integer(i - 1L), 0L),
            FUN.VALUE = character(1)
        )
    )
})

test_that("as.matrix() rejects an unknown reading", {
    expect_error(as.matrix(tasGeno, type = "probability"), "should be one of")
})

test_that("a genomic dataset passes 'type' to its genotype half", {
    expect_equal(
        as.matrix(tasDataset[, sites(1:10)], type = "allele"),
        as.matrix(genotype(tasDataset)[, sites(1:10)], type = "allele")
    )
})


## Numeric genotype matrix ----

test_that("as.matrix() on a numeric genotype returns reference probabilities", {
    numGt <- readGenotype(returnSysFiles("numeric_genotype.txt"))
    m <- as.matrix(numGt)

    expect_type(m, "double")
    expect_equal(dim(m), c(3, 3))
    expect_equal(rownames(m), taxaList(numGt))
    expect_equal(colnames(m), positionList(numGt)$Name)

    jGt <- javaRefObj(numGt)
    expect_equal(m[1, 1], jGt$referenceProbability(0L, 0L))
    expect_equal(m[3, 2], jGt$referenceProbability(2L, 1L))
})

test_that("a numeric genotype read from an R matrix round trips", {
    m <- rtMatrices$num_gt_md
    back <- as.matrix(readGenotype(m))

    expect_equal(dimnames(back), dimnames(m))

    # TASSEL stores a reference probability in one byte, so the values
    # come back quantized to 1/255 of the range rather than unchanged
    expect_lt(max(abs(back - m)), 1 / 255)
})


## SummarizedExperiment coercion ----

test_that("as(<genotype>, 'SummarizedExperiment') puts sites in rows", {
    sub <- tasGeno[, sites(1:25)]
    se <- methods::as(sub, "SummarizedExperiment")

    expect_s4_class(se, "RangedSummarizedExperiment")
    expect_equal(dim(se), c(25, 281))
    expect_equal(
        unname(SummarizedExperiment::assay(se)),
        unname(t(as.matrix(sub)))
    )
    expect_equal(rownames(se), positionList(sub)$Name)
    expect_equal(colnames(se), taxaList(sub))
    expect_equal(
        SummarizedExperiment::rowRanges(se),
        granges(sub, use.mcols = TRUE)
    )
})

test_that("as(<dataset>, 'SummarizedExperiment') takes the genotype half", {
    sub <- tasDataset[, sites(1:25)]
    se <- methods::as(sub, "SummarizedExperiment")

    expect_equal(dim(se), c(25, length(taxaList(sub))))
    expect_equal(colnames(se), taxaList(sub))
})


## granges() ----

test_that("granges() returns one width-1 range per marker", {
    gr <- granges(tasGeno)
    pl <- positionList(tasGeno)

    expect_s4_class(gr, "GRanges")
    expect_equal(length(gr), nrow(pl))
    expect_true(all(GenomicRanges::width(gr) == 1L))
    expect_equal(names(gr), pl$Name)
    expect_equal(as.character(GenomicRanges::seqnames(gr)), pl$Chromosome)
    expect_equal(GenomicRanges::start(gr), pl$Position)
})

test_that("granges() agrees with seqnames() on the chromosomes present", {
    expect_setequal(
        levels(GenomicRanges::seqnames(granges(tasGeno))),
        seqnames(tasGeno)
    )
})

test_that("granges() carries the remaining position columns on request", {
    bare  <- granges(tasGeno)
    mcold <- granges(tasGeno, use.mcols = TRUE)

    expect_equal(ncol(S4Vectors::mcols(bare)), 0L)
    expect_equal(
        colnames(S4Vectors::mcols(mcold)),
        c("Site", "Name", "VARIANT")
    )
    expect_equal(S4Vectors::mcols(mcold)$Site, positionList(tasGeno)$Site)
})

test_that("granges() can drop the marker names", {
    expect_null(names(granges(tasGeno, use.names = FALSE)))
})

test_that("granges() is unaffected by the class holding the positions", {
    expect_equal(granges(tasGeno), granges(tasDataset))
    expect_equal(granges(tasGeno), granges(rtObjsLegacy$gt_hmp))
})

test_that("granges() output can be filtered on", {
    target <- granges(tasGeno[, region("2", 20e6, 40e6)])
    sub <- tasGeno[, sitesWhere(overlaps(target))]

    expect_equal(nrow(positionList(sub)), length(target))
    expect_equal(positionList(sub)$Name, names(target))
})


## Distance matrix ----

test_that("as.matrix() on a distance matrix labels both axes", {
    kin <- kinshipMatrix(tasGeno)
    m <- as.matrix(kin)

    expect_true(is.matrix(m))
    expect_type(m, "double")
    expect_equal(dim(m), c(281, 281))
    expect_equal(rownames(m), taxaList(kin))
    expect_equal(colnames(m), taxaList(kin))
})

test_that("the array route agrees with the tab-delimited text it replaced", {
    kin <- kinshipMatrix(tasGeno)

    # How 'as.matrix()' read a distance matrix before 0.14.0. The text
    # form rounds, so the two agree only to the precision it carries
    parsed <- local({
        rows <- unlist(strsplit(javaRefObj(kin)$toStringTabDelim(), split = "\n"))
        cells <- t(simplify2array(strsplit(rows, split = "\t")))
        colnames(cells) <- as.character(unlist(cells[1, ]))
        cells <- cells[-1, ]
        taxa <- cells[, 1]
        out <- apply(cells[, -1], 2, as.numeric)
        rownames(out) <- taxa
        out
    })

    expect_equal(as.matrix(kin), parsed, tolerance = 1e-6)
    expect_equal(dimnames(as.matrix(kin)), dimnames(parsed))
})

test_that("a distance matrix answers javaRefObj() and taxaList()", {
    kin <- kinshipMatrix(tasGeno)

    expect_s4_class(javaRefObj(kin), "jobjRef")
    expect_true(
        javaRefObj(kin) %instanceof% "net.maizegenetics.taxa.distance.DistanceMatrix"
    )
    expect_equal(taxaList(kin), rownames(kin))
    expect_length(taxaList(kin), 281)
})

test_that("as.dist() returns the lower triangle of the pairwise matrix", {
    kin <- kinshipMatrix(tasGeno)
    d <- as.dist(kin)

    expect_s3_class(d, "dist")
    expect_equal(attr(d, "Size"), 281L)
    expect_equal(labels(d), taxaList(kin))

    # A 'dist' keeps only the lower triangle, so the self-similarity on
    # the diagonal of a kinship matrix is dropped
    square <- as.matrix(kin)
    expect_equal(as.matrix(d)[lower.tri(square)], square[lower.tri(square)])
    expect_true(all(diag(as.matrix(d)) == 0))
})

test_that("as.dist() output is accepted by hclust()", {
    dm <- distanceMatrix(tasGeno[, sites(1:200)])

    expect_s3_class(stats::hclust(as.dist(dm)), "hclust")
})


## Java string helper ----

test_that(".jStrings() resolves a list of Java objects", {
    chroms <- javaRefObj(tasGeno)$chromosomes()

    expect_equal(.jStrings(chroms), seqnames(tasGeno))
    expect_type(.jStrings(chroms), "character")
    expect_equal(.jStrings(list()), character(0))
})


## Table report conventions ----

test_that("tableReportToDF() replaces spaces in column names", {
    siteDf <- siteSummary(tasGeno)

    expect_false(any(grepl(" ", colnames(siteDf))))
    expect_true("Minor_Allele_Frequency" %in% colnames(siteDf))
})

test_that("tableReportToDF() spells missing values by column type", {
    phenoDf <- as.data.frame(tasPheno)

    # A missing number arrives as NaN, which 'is.na()' still recognises
    expect_true(any(is.nan(phenoDf$EarHT)))
    expect_true(all(is.na(phenoDf$EarHT[is.nan(phenoDf$EarHT)])))
})


## Attribute metadata schema ----

test_that("the trait metadata helpers agree on their column names", {
    internal <- extractPhenotypeAttDf(javaRefObj(tasPheno))
    public <- attributeData(tasPheno)

    expect_equal(
        colnames(internal),
        c("trait_id", "trait_type", "trait_attribute")
    )
    expect_equal(internal$trait_id, public$trait_id)
    expect_equal(internal$trait_type, public$trait_type)
    expect_equal(internal$trait_attribute, public$trait_attribute)
})


## Phenotype round trip ----

test_that("attributeData() can be handed back to readPhenotype()", {
    back <- readPhenotype(
        as.data.frame(tasPheno),
        attr = attributeData(tasPheno)
    )

    expect_s4_class(back, "TasselPhenotype")
    expect_equal(as.data.frame(back), as.data.frame(tasPheno))
    expect_equal(traitNames(back), traitNames(tasPheno))
    expect_equal(
        attributeData(back)$trait_type,
        attributeData(tasPheno)$trait_type
    )
})

test_that("readPhenotype() still takes the col_id spelling", {
    attrDf <- tibble::tribble(
        ~"col_id",      ~"tassel_attr",
        "taxa_id",      "taxa",
        "plant_height", "data",
        "PC1",          "covariate"
    )
    df <- tibble::tribble(
        ~"taxa_id", ~"plant_height", ~"PC1",
        "line_a",   12.3,            0.5,
        "line_b",   22.8,            -1.5
    )

    out <- readPhenotype(df, attr = attrDf)

    expect_s4_class(out, "TasselPhenotype")
    expect_equal(traitNames(out), c("plant_height", "PC1"))
})

test_that("validateAttrDf() normalises either spelling and rejects neither", {
    attrData <- tibble::tibble(trait_id = "a", trait_type = "taxa")
    colId    <- tibble::tibble(col_id = "a", tassel_attr = "taxa")

    expect_equal(names(validateAttrDf(attrData)), names(colId))
    expect_equal(validateAttrDf(colId), colId)
    expect_error(validateAttrDf(data.frame(a = 1)), "Incorrect column IDs")
    expect_error(validateAttrDf(list()), "needs to be of type 'data.frame'")
})


## tableReport() contract ----

test_that("tableReport() reads a report name the same way on every class", {
    ldRes <- linkageDiseq(
        tasGeno[, region("2", 228e6, 300e6)],
        ldType  = "All",
        verbose = FALSE
    )

    expect_equal(reportNames(ldRes), "LD")
    expect_s3_class(tableReport(ldRes), "tbl_df")
    expect_equal(tableReport(ldRes, "LD"), tableReport(ldRes))
    expect_equal(names(tableReport(ldRes, "ALL")), "LD")
    expect_error(tableReport(ldRes, "NOPE"), "Report ID not found")
})

test_that("the base AssociationResults class accepts the catch all", {
    obj <- methods::new(
        "AssociationResults",
        results   = list(one = data.frame(x = 1), two = data.frame(y = 2)),
        traits    = "trait_1",
        assocType = "GLM"
    )

    # No single default report, so a missing name means every report
    expect_equal(tableReport(obj), tableReport(obj, "ALL"))
    expect_equal(names(tableReport(obj, "ALL")), c("one", "two"))
    expect_equal(tableReport(obj, "one"), data.frame(x = 1))
    expect_error(tableReport(obj, "NOPE"), "Report ID not found")
    expect_error(tableReport(obj, 1), "must be of type 'character'")
})
