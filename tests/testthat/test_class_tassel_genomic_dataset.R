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

    dsDf <- readGenomicDataset(gt, phDf, attrTypes = attrDf)

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

test_that("readGenomicDataset() requires attrTypes for a data frame phenotype", {
    expect_error(
        readGenomicDataset(gt, data.frame(taxa = "33-16", yield = 1)),
        "needs attribute metadata"
    )
})

test_that("readGenomicDataset() reports a join with no shared taxa", {
    phUnrelated <- readPhenotype(
        data.frame(taxa_id = c("zz1", "zz2"), yield = c(1, 2)),
        attrTypes = data.frame(
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

test_that("show() reports dimensions, columns, and the Java address", {
    out <- capture_output(show(ds))

    expect_match(out, "TasselGenomicDataset")
    expect_match(out, "278 taxa")
    expect_match(out, "3093 sites")
    expect_match(out, "Column types: taxa: 1, data: 3, genotype: 1")
    expect_match(out, ds@jMemAddress)
})

test_that("show() prints the phenotype traits as tibble columns", {
    out <- capture_output(show(ds))

    expect_match(out, "Taxa")
    expect_match(out, "EarHT")
    expect_match(out, "<data>")
    expect_match(out, "33-16")
})

test_that("show() adds a Genotype column holding the leading sites", {
    out <- capture_output(show(ds))

    expect_match(out, "Genotype")
    expect_match(out, "<geno>")
    expect_match(out, "Genotype: showing the first 5 of 3093 sites")
})

test_that("show() reports the number of observations when truncated", {
    out <- capture_output(show(rtObjs$ds_hmp_ph_full))

    expect_match(out, "showing the first 10 of 525 rows")
})

test_that("show() breaks columns down by type", {
    out <- capture_output(show(rtObjs$ds_hmp_ph_full))

    expect_match(
        out,
        "Column types: taxa: 1, factor: 1, data: 3, covariate: 3, genotype: 1"
    )
})

test_that("show() keeps the Genotype column on a narrow console", {
    dispData <- rTASSEL:::formatGenomicDatasetDisplay(rtObjs$ds_hmp_ph_full)
    out      <- paste(capture_output(print(dispData, width = 40)), collapse = "\n")

    # The genotype column holds its place while the trailing traits are
    # squeezed into pillar's own footer
    expect_match(out, "<geno>")
    expect_match(out, "more variables")
    expect_match(out, "Q3 <cov>")
})

test_that("show() omits the site notice when every site is displayed", {
    out <- capture_output(show(ds[, sites(1:3)]))

    expect_match(out, "3 sites")
    expect_false(grepl("showing the first 3 of 3 sites", out))
})

test_that("show() displays reference probabilities for numeric genotypes", {
    dsNum <- readGenomicDataset(
        rtMatrices$num_gt_sm,
        data.frame(Taxa = rownames(rtMatrices$num_gt_sm), yield = c(1, 2, 3)),
        attrTypes = data.frame(
            col_id      = c("Taxa", "yield"),
            tassel_attr = c("taxa", "data")
        )
    )

    out <- capture_output(show(dsNum))

    expect_match(out, "<geno>")
    expect_match(out, "0\\.\\d{3}")
    expect_match(out, "Column types: taxa: 1, data: 1, genotype: 1")
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
