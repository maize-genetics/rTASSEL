# === Tests for the TasselGenomicDataset class =======================

## Preamble ----
startLogger()

ds       <- rtObjs$ds_hmp_ph_nomiss
gt       <- rtObjs$gt_hmp
phNoMiss <- rtObjs$ph_nomiss

# The hapmap file holds 281 taxa and the "no missing" phenotype file 298;
# 278 taxa are shared between them
nGtTaxa  <- 281L
nPhTaxa  <- 298L
nBothTaxa <- 278L
nSites   <- 3093L


# /// readGenomicDataset() ///////////////////////////////////////////

test_that("readGenomicDataset() builds a dataset from file paths", {
    dsPaths <- readGenomicDataset(rtFiles$gt_hmp_path, rtFiles$ph_nomiss_path)

    expect_s4_class(dsPaths, "TasselGenomicDataset")
    expect_equal(dsPaths@genotype@jRefObj$numberOfTaxa(), nBothTaxa)
    expect_equal(dsPaths@genotype@jRefObj$numberOfSites(), nSites)
})

test_that("readGenomicDataset() builds a dataset from rTASSEL objects", {
    expect_s4_class(ds, "TasselGenomicDataset")
    expect_equal(ds@jClass, "net.maizegenetics.phenotype.GenotypePhenotype")
    expect_true(nzchar(ds@jMemAddress))
})

test_that("readGenomicDataset() builds a dataset from a data frame", {
    attrDf <- data.frame(
        col_id      = c("Taxa", "yield"),
        tassel_attr = c("taxa", "data")
    )
    phDf <- data.frame(
        Taxa  = head(taxaList(gt), 5),
        yield = c(1, 2, 3, 4, 5)
    )

    dsDf <- readGenomicDataset(gt, phDf, attr = attrDf)

    expect_s4_class(dsDf, "TasselGenomicDataset")
    expect_equal(dsDf@genotype@jRefObj$numberOfTaxa(), 5L)
    expect_equal(traitNames(dsDf), "yield")
})

test_that("readGenomicDataset() mixes an object with a path", {
    dsMixed <- readGenomicDataset(gt, rtFiles$ph_nomiss_path)

    expect_s4_class(dsMixed, "TasselGenomicDataset")
    expect_equal(dsMixed@genotype@jRefObj$numberOfTaxa(), nBothTaxa)
})


# /// Join semantics /////////////////////////////////////////////////

test_that("readGenomicDataset() intersects taxa by default", {
    expect_equal(ds@genotype@jRefObj$numberOfTaxa(), nBothTaxa)
    expect_equal(nrow(as.data.frame(ds)), nBothTaxa)
})

test_that("readGenomicDataset() keeps all taxa when join = 'union'", {
    dsUnion <- readGenomicDataset(gt, phNoMiss, join = "union")

    expect_equal(dsUnion@genotype@jRefObj$numberOfTaxa(), nGtTaxa)
    expect_equal(nrow(as.data.frame(dsUnion)), nPhTaxa)
})

test_that("readGenomicDataset() rejects an unknown join type", {
    expect_error(readGenomicDataset(gt, phNoMiss, join = "outer"), "must be one of")
})


# /// Constructor errors /////////////////////////////////////////////

test_that("readGenomicDataset() rejects inputs missing genotype data", {
    expect_error(
        readGenomicDataset(phNoMiss, phNoMiss),
        "`genotype` does not contain genotype data"
    )
})

test_that("readGenomicDataset() rejects inputs missing phenotype data", {
    expect_error(
        readGenomicDataset(gt, gt),
        "`phenotype` does not contain phenotype data"
    )
})

test_that("readGenomicDataset() requires attr for a data frame phenotype", {
    expect_error(
        readGenomicDataset(gt, data.frame(taxa = "33-16", yield = 1)),
        "need attribute metadata"
    )
})

test_that("readGenomicDataset() reports a join with no shared taxa", {
    phUnrelated <- readPhenotype(
        data.frame(taxa_id = c("zz1", "zz2"), yield = c(1, 2)),
        attr = data.frame(
            col_id      = c("taxa_id", "yield"),
            tassel_attr = c("taxa", "data")
        )
    )

    expect_error(
        readGenomicDataset(gt, phUnrelated),
        "Could not join the genotype and phenotype data"
    )
})

test_that("createTasselGenomicDataset() validates its Java input", {
    expect_error(
        createTasselGenomicDataset(rJava::.jnull()),
        "must be a non-null Java object reference"
    )
    expect_error(
        createTasselGenomicDataset(getGenotypeTable(gt)),
        "must hold both genotype and phenotype data"
    )
})


# /// Component accessors ////////////////////////////////////////////

test_that("genotype() and phenotype() return the joined components", {
    expect_s4_class(genotype(ds), "TasselGenotype")
    expect_s4_class(phenotype(ds), "TasselPhenotype")

    # The components wrap the joined tables, not the originals handed in
    expect_equal(genotype(ds)@jRefObj$numberOfTaxa(), nBothTaxa)
    expect_lt(genotype(ds)@jRefObj$numberOfTaxa(), gt@jRefObj$numberOfTaxa())
})

test_that("javaRefObj() returns the backing GenotypePhenotype", {
    expect_equal(
        rJava::.jclass(javaRefObj(ds)),
        "net.maizegenetics.phenotype.GenotypePhenotype"
    )
})


# /// Taxa and position methods //////////////////////////////////////

test_that("taxaList() reflects the joined taxa", {
    expect_type(taxaList(ds), "character")
    expect_length(taxaList(ds), nBothTaxa)
})

test_that("positionList() matches the genotype component", {
    expect_equal(nrow(positionList(ds)), nSites)
    expect_equal(positionList(ds), positionList(genotype(ds)))
})

test_that("seqnames() returns the dataset's chromosomes", {
    expect_equal(seqnames(ds), as.character(1:10))
})


# /// Summary methods ////////////////////////////////////////////////

test_that("siteSummary() and taxaSummary() delegate to the genotype", {
    expect_equal(nrow(siteSummary(ds)), nSites)
    expect_equal(nrow(taxaSummary(ds)), nBothTaxa)
})

test_that("attributeData() and traitNames() delegate to the phenotype", {
    expect_equal(attributeData(ds), attributeData(phenotype(ds)))
    expect_equal(traitNames(ds), c("EarHT", "dpoll", "EarDia"))
})

test_that("attributeData() carries every phenotype attribute type", {
    attrDf <- attributeData(rtObjs$ds_hmp_ph_full)

    expect_equal(
        attrDf$trait_type,
        c("taxa", "factor", "data", "data", "data", rep("covariate", 3))
    )
})


# /// show ///////////////////////////////////////////////////////////

test_that("show() reports dimensions, traits, and the Java address", {
    out <- capture_output(show(ds))

    expect_match(out, "TasselGenomicDataset")
    expect_match(out, "278 taxa")
    expect_match(out, "3093 sites")
    expect_match(out, "<TasselGenotype>")
    expect_match(out, "4 traits \\(data: 3, taxa: 1\\)")
    expect_match(out, ds@jMemAddress)
})

test_that("show() breaks traits down by attribute type", {
    out <- capture_output(show(rtObjs$ds_hmp_ph_full))

    expect_match(out, "8 traits \\(covariate: 3, data: 3, factor: 1, taxa: 1\\)")
})

test_that("formatTraitSummary() handles a phenotype with no traits", {
    expect_equal(formatTraitSummary(list()), "no traits")
    expect_equal(formatTraitSummary(list(data = 2, taxa = 1)), "data: 2, taxa: 1")
})


# /// Coercion ///////////////////////////////////////////////////////

test_that("as.data.frame() returns the joined phenotype data", {
    df <- as.data.frame(ds)

    expect_s3_class(df, "tbl_df")
    expect_equal(dim(df), c(nBothTaxa, 4L))
    expect_equal(colnames(df), c("Taxa", "EarHT", "dpoll", "EarDia"))
    expect_equal(df, as.data.frame(phenotype(ds)))
})

test_that("as.matrix() returns the joined genotype dosages", {
    m <- as.matrix(ds)

    expect_true(is.matrix(m))
    expect_equal(dim(m), c(nBothTaxa, nSites))
    expect_equal(rownames(m), taxaList(ds))
})


# /// Bracket subsetting /////////////////////////////////////////////

test_that("[ keeps both components in step when selecting taxa", {
    sub <- ds[taxa("33-16", "38-11"), ]

    expect_s4_class(sub, "TasselGenomicDataset")
    expect_equal(sub@genotype@jRefObj$numberOfTaxa(), 2L)
    expect_equal(nrow(as.data.frame(sub)), 2L)
    expect_setequal(as.data.frame(sub)$Taxa, c("33-16", "38-11"))
})

test_that("[ leaves taxa untouched when only sites are selected", {
    sub <- ds[, siteIds("PZB00859.1", "PZA01271.1")]

    expect_equal(sub@genotype@jRefObj$numberOfSites(), 2L)
    expect_equal(sub@genotype@jRefObj$numberOfTaxa(), nBothTaxa)
    expect_equal(nrow(as.data.frame(sub)), nBothTaxa)
})

test_that("[ applies taxa and site selectors together", {
    sub <- ds[taxa("33-16", "38-11", "4226"), chrom("5")]

    expect_equal(sub@genotype@jRefObj$numberOfTaxa(), 3L)
    expect_equal(
        sub@genotype@jRefObj$numberOfSites(),
        gt[, chrom("5")]@jRefObj$numberOfSites()
    )
})

test_that("[ without selectors returns an equivalent dataset", {
    sub <- ds[]

    expect_s4_class(sub, "TasselGenomicDataset")
    expect_equal(sub@genotype@jRefObj$numberOfTaxa(), nBothTaxa)
    expect_equal(sub@genotype@jRefObj$numberOfSites(), nSites)
})


# /// Downstream use /////////////////////////////////////////////////

test_that("a dataset can drive an analysis that needs both data types", {
    kin <- kinshipMatrix(ds)

    expect_s4_class(kin, "TasselDistanceMatrix")
    expect_equal(dim(as.matrix(kin)), c(nBothTaxa, nBothTaxa))
})

test_that(".wrapLikeInput() rebuilds a dataset from a Java result", {
    rebuilt <- rTASSEL:::.wrapLikeInput(javaRefObj(ds), ds)

    expect_s4_class(rebuilt, "TasselGenomicDataset")
    expect_equal(rebuilt@genotype@jRefObj$numberOfTaxa(), nBothTaxa)
})
