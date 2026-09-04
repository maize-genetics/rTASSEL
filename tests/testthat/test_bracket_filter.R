# === Tests for bracket-based filtering ==============================

## Preamble ----
startLogger()

gt <- readGenotype(rtFiles$gt_hmp_path)
ds <- readGenomicDataset(gt, rtFiles$ph_nomiss_path)

testGr <- GenomicRanges::GRanges(
    seqnames = c("1", "2"),
    ranges   = IRanges::IRanges(
        start = c(1e6, 5e5),
        end   = c(5e6, 1e7)
    )
)


# /// Selector construction tests ////////////////////////////////////

test_that("taxa() creates a valid TaxaSelector", {
    sel <- taxa("33-16", "38-11")
    expect_s4_class(sel, "TaxaSelector")
    expect_equal(sel@type, "ids")
    expect_equal(sel@ids, c("33-16", "38-11"))
    expect_false(sel@negate)
})

test_that("taxaWhere() creates a valid predicate TaxaSelector", {
    sel <- taxaWhere(startsWith(taxaId, "A"))
    expect_s4_class(sel, "TaxaSelector")
    expect_equal(sel@type, "predicate")
    expect_false(sel@negate)
})

test_that("sites() creates a valid SiteSelector", {
    sel <- sites(1:100)
    expect_s4_class(sel, "SiteSelector")
    expect_equal(sel@type, "indices")
    expect_equal(sel@indices, 1:100L)
    expect_false(sel@negate)
})

test_that("sites() rejects indices below 1", {
    expect_error(sites(0:99), "1-based")
    expect_error(sites(-1), "1-based")
})

test_that("siteIds() creates a valid SiteSelector", {
    sel <- siteIds("PZB00859.1", "PZA01271.1")
    expect_s4_class(sel, "SiteSelector")
    expect_equal(sel@type, "names")
    expect_equal(sel@ids, c("PZB00859.1", "PZA01271.1"))
})

test_that("chrom() creates a valid SiteSelector", {
    sel <- chrom("5")
    expect_s4_class(sel, "SiteSelector")
    expect_equal(sel@type, "chrom")
    expect_equal(sel@chromId, "5")
})

test_that("region() creates a valid SiteSelector", {
    sel <- region("1", 1e6, 2e6)
    expect_s4_class(sel, "SiteSelector")
    expect_equal(sel@type, "region")
    expect_equal(sel@chromId, "1")
    expect_equal(sel@start, 1e6)
    expect_equal(sel@end, 2e6)
})

test_that("region() creates a granges SiteSelector from a GRanges", {
    sel <- region(testGr)
    expect_s4_class(sel, "SiteSelector")
    expect_equal(sel@type, "granges")
    expect_identical(sel@granges, testGr)
})

test_that("region() rejects coordinates alongside a GRanges", {
    expect_error(region(testGr, 1, 100))
    expect_error(region(GenomicRanges::GRanges()))
})

test_that("sitesWhere() creates a valid predicate SiteSelector", {
    sel <- sitesWhere(maf >= 0.05)
    expect_s4_class(sel, "SiteSelector")
    expect_equal(sel@type, "predicate")
    expect_false(sel@negate)
})


# /// Negation tests /////////////////////////////////////////////////

test_that("! toggles negate on TaxaSelector", {
    sel <- taxa("33-16")
    neg <- !sel
    expect_true(neg@negate)
    expect_false((!neg)@negate)
})

test_that("! toggles negate on SiteSelector", {
    sel <- sites(1:10)
    neg <- !sel
    expect_true(neg@negate)
    expect_false((!neg)@negate)
})


# /// Constructor error tests ////////////////////////////////////////

test_that("constructors reject empty input", {
    expect_error(taxa())
    expect_error(sites())
    expect_error(siteIds())
    expect_error(chrom())
})


# /// Bracket filtering: taxa ////////////////////////////////////////

test_that("[taxa()] selects specific taxa by ID", {
    ids <- c("33-16", "38-11", "4226")
    sub <- gt[taxa(ids), ]
    expect_s4_class(sub, "TasselGenotype")
    expect_equal(sub@jRefObj$numberOfTaxa(), length(ids))
})

test_that("[character] selects taxa by plain character vector", {
    ids <- c("33-16", "38-11")
    sub <- gt[ids, ]
    expect_s4_class(sub, "TasselGenotype")
    expect_equal(sub@jRefObj$numberOfTaxa(), 2L)
})

test_that("[taxaWhere()] filters taxa by predicate", {
    sub <- gt[taxaWhere(startsWith(taxaId, "A")), ]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfTaxa() > 0)
    expect_true(sub@jRefObj$numberOfTaxa() < gt@jRefObj$numberOfTaxa())
})

test_that("[!taxa()] negates taxa selection", {
    ids <- c("33-16", "38-11")
    sub <- gt[!taxa(ids), ]
    totalTaxa <- gt@jRefObj$numberOfTaxa()
    expect_equal(sub@jRefObj$numberOfTaxa(), totalTaxa - length(ids))
})


# /// Taxa metrics ///////////////////////////////////////////////////

test_that("buildTaxaMetadata() returns taxa IDs and proportions", {
    meta <- rTASSEL:::buildTaxaMetadata(gt@jRefObj)
    expect_named(meta, c("taxaId", "notMissing", "het"))
    expect_equal(nrow(meta), gt@jRefObj$numberOfTaxa())
    expect_true(all(meta$notMissing >= 0 & meta$notMissing <= 1))
    expect_true(all(meta$het >= 0 & meta$het <= 1))
})

test_that("buildTaxaMetadata() computes only the requested columns", {
    meta <- rTASSEL:::buildTaxaMetadata(gt@jRefObj, needed = "taxaId")
    expect_named(meta, "taxaId")
})

test_that("[taxaWhere()] filters taxa by notMissing", {
    sub <- gt[taxaWhere(notMissing >= 0.9), ]
    meta <- rTASSEL:::buildTaxaMetadata(gt@jRefObj, needed = "notMissing")
    expect_equal(sub@jRefObj$numberOfTaxa(), sum(meta$notMissing >= 0.9))
})

test_that("[taxaWhere()] filters taxa by het", {
    sub <- gt[taxaWhere(het <= 0.01), ]
    meta <- rTASSEL:::buildTaxaMetadata(gt@jRefObj, needed = "het")
    expect_equal(sub@jRefObj$numberOfTaxa(), sum(meta$het <= 0.01))
})

test_that("taxa plugin short-circuit agrees with the R fallback", {
    meta <- rTASSEL:::buildTaxaMetadata(gt@jRefObj)
    expected <- sum(meta$notMissing >= 0.9 & meta$het <= 0.01)

    # The literal form is mapped onto FilterTaxaBuilderPlugin, the form
    # referencing an outer variable is not, so both paths get exercised
    thresh <- 0.9
    shortCircuited <- gt[taxaWhere(notMissing >= 0.9 & het <= 0.01), ]
    fellBack       <- gt[taxaWhere(notMissing >= thresh & het <= 0.01), ]

    expect_equal(shortCircuited@jRefObj$numberOfTaxa(), expected)
    expect_equal(fellBack@jRefObj$numberOfTaxa(), expected)
})

test_that("taxaWhere() metrics match the legacy taxa filter", {
    withr::local_options(lifecycle_verbosity = "quiet")
    legacy <- filterGenotypeTableTaxa(
        readGenotypeTableFromPath(rtFiles$gt_hmp_path),
        minNotMissing = 0.9
    )
    sub <- gt[taxaWhere(notMissing >= 0.9), ]
    expect_equal(sub@jRefObj$numberOfTaxa(), getTaxaList(legacy)$size())
})

test_that("[!taxaWhere()] negates a metric predicate", {
    meta <- rTASSEL:::buildTaxaMetadata(gt@jRefObj, needed = "het")
    sub <- gt[!taxaWhere(het <= 0.01), ]
    expect_equal(sub@jRefObj$numberOfTaxa(), sum(meta$het > 0.01))
})


# /// Bracket filtering: sites ///////////////////////////////////////

test_that("[sites()] selects by 1-based index", {
    sub <- gt[, sites(1:100)]
    expect_s4_class(sub, "TasselGenotype")
    expect_equal(sub@jRefObj$numberOfSites(), 100L)
})

test_that("[sites()] indices line up with positionList()", {
    pl <- positionList(gt)
    sub <- gt[, sites(c(1, 7, 42))]
    expect_equal(
        rTASSEL:::batchSiteNames(sub@jRefObj),
        pl$Name[c(1, 7, 42)]
    )
})

test_that("[numeric] selects sites by plain numeric vector", {
    sub <- gt[, 1:50]
    expect_s4_class(sub, "TasselGenotype")
    expect_equal(sub@jRefObj$numberOfSites(), 50L)
})

test_that("[numeric] rejects indices below 1", {
    expect_error(gt[, 0:49], "1-based")
})

test_that("[sitesWhere()] exposes 1-based siteIndex", {
    sub <- gt[, sitesWhere(siteIndex <= 10)]
    expect_equal(sub@jRefObj$numberOfSites(), 10L)
    expect_equal(
        rTASSEL:::batchSiteNames(sub@jRefObj),
        positionList(gt)$Name[1:10]
    )
})

test_that("[siteIds()] selects by marker name", {
    names <- c("PZB00859.1", "PZA01271.1")
    sub <- gt[, siteIds(names)]
    expect_s4_class(sub, "TasselGenotype")
    expect_equal(sub@jRefObj$numberOfSites(), length(names))
})

test_that("[chrom()] selects all sites on a chromosome", {
    sub <- gt[, chrom("5")]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfSites() > 0)
    expect_true(sub@jRefObj$numberOfSites() < gt@jRefObj$numberOfSites())
})

test_that("[chrom()] selects multiple chromosomes", {
    sub <- gt[, chrom("5", "10")]
    subSingle5  <- gt[, chrom("5")]
    subSingle10 <- gt[, chrom("10")]
    expect_equal(
        sub@jRefObj$numberOfSites(),
        subSingle5@jRefObj$numberOfSites() + subSingle10@jRefObj$numberOfSites()
    )
})

test_that("[region()] selects a genomic region", {
    sub <- gt[, region("1", 1e6, 5e6)]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfSites() > 0)
    expect_true(sub@jRefObj$numberOfSites() < gt@jRefObj$numberOfSites())
})

test_that("[region()] selects sites from a GRanges object", {
    sub <- gt[, region(testGr)]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfSites() > 0)
    expect_true(sub@jRefObj$numberOfSites() < gt@jRefObj$numberOfSites())
})

test_that("[region(GRanges)] matches the legacy gRangesObj filter", {
    withr::local_options(lifecycle_verbosity = "quiet")
    legacy <- filterGenotypeTableSites(
        readGenotypeTableFromPath(rtFiles$gt_hmp_path),
        gRangesObj = testGr
    )
    sub <- gt[, region(testGr)]
    expect_equal(
        sub@jRefObj$numberOfSites(),
        getGenotypeTable(legacy)$numberOfSites()
    )
})

test_that("[!region(GRanges)] negates a GRanges selection", {
    nTotal <- gt@jRefObj$numberOfSites()
    nKept  <- gt[, region(testGr)]@jRefObj$numberOfSites()
    sub <- gt[, !region(testGr)]
    expect_equal(sub@jRefObj$numberOfSites(), nTotal - nKept)
})

test_that("[sitesWhere()] filters sites by MAF predicate", {
    sub <- gt[, sitesWhere(maf >= 0.05)]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfSites() > 0)
    expect_true(sub@jRefObj$numberOfSites() < gt@jRefObj$numberOfSites())
})

test_that("[sitesWhere()] filters by compound predicate", {
    sub <- gt[, sitesWhere(chrom == "1" & maf >= 0.05)]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfSites() > 0)

    subChromOnly <- gt[, chrom("1")]
    expect_true(sub@jRefObj$numberOfSites() <= subChromOnly@jRefObj$numberOfSites())
})

test_that("[sitesWhere()] filters sites by alleleCount predicate", {
    sub <- gt[, sitesWhere(alleleCount >= 10)]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfSites() > 0)
    expect_true(sub@jRefObj$numberOfSites() <= gt@jRefObj$numberOfSites())
})

test_that("[sitesWhere()] filters sites by het predicate", {
    sub <- gt[, sitesWhere(het <= 0.5)]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfSites() > 0)
    expect_true(sub@jRefObj$numberOfSites() <= gt@jRefObj$numberOfSites())
})

test_that("[sitesWhere()] filters sites by isIndel predicate", {
    sub <- gt[, sitesWhere(!isIndel)]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfSites() > 0)
    expect_true(sub@jRefObj$numberOfSites() <= gt@jRefObj$numberOfSites())
})

test_that("[sitesWhere()] filters sites by isBiallelic predicate", {
    sub <- gt[, sitesWhere(isBiallelic)]
    expect_s4_class(sub, "TasselGenotype")
    expect_true(sub@jRefObj$numberOfSites() > 0)
    expect_true(sub@jRefObj$numberOfSites() <= gt@jRefObj$numberOfSites())
})

test_that("[!sites()] negates site selection", {
    nTotal <- gt@jRefObj$numberOfSites()
    sub <- gt[, !sites(1:10)]
    expect_equal(sub@jRefObj$numberOfSites(), nTotal - 10L)
})


# /// Combined taxa + sites //////////////////////////////////////////

test_that("bracket filters taxa and sites simultaneously", {
    sub <- gt[taxa("33-16", "38-11"), sites(1:50)]
    expect_equal(sub@jRefObj$numberOfTaxa(), 2L)
    expect_equal(sub@jRefObj$numberOfSites(), 50L)
})

test_that("combined taxa + sitesWhere predicate", {
    sub <- gt[taxa("33-16"), sitesWhere(maf >= 0.05)]
    expect_equal(sub@jRefObj$numberOfTaxa(), 1L)
    expect_true(sub@jRefObj$numberOfSites() > 0)
})


# /// Return type preservation ///////////////////////////////////////

test_that("bracket returns same class as input (TasselGenotype)", {
    sub <- gt[taxa("33-16"), ]
    expect_s4_class(sub, "TasselGenotype")
    expect_false(is(sub, "TasselNumericGenotype"))
})

test_that("bracket returns TasselNumericGenotype when input is", {
    numGt <- readGenotype(rtMatrices$num_gt_md)
    expect_s4_class(numGt, "TasselNumericGenotype")
    sub <- numGt[, sites(1:5)]
    expect_s4_class(sub, "TasselNumericGenotype")
})


# /// Bracket filtering on a genomic dataset /////////////////////////

test_that("bracket subsets a TasselGenomicDataset by sites", {
    sub <- ds[, sites(1:50)]
    expect_s4_class(sub, "TasselGenomicDataset")
    expect_equal(sub@genotype@jRefObj$numberOfSites(), 50L)
    expect_equal(
        sub@genotype@jRefObj$numberOfTaxa(),
        ds@genotype@jRefObj$numberOfTaxa()
    )
})

test_that("bracket subsets a TasselGenomicDataset by a GRanges region", {
    sub <- ds[, region(testGr)]
    expect_s4_class(sub, "TasselGenomicDataset")
    expect_equal(
        sub@genotype@jRefObj$numberOfSites(),
        gt[, region(testGr)]@jRefObj$numberOfSites()
    )
})

test_that("bracket re-joins phenotype data after a taxa predicate", {
    sub <- ds[taxaWhere(notMissing >= 0.9), ]
    expect_s4_class(sub, "TasselGenomicDataset")
    expect_true(
        sub@genotype@jRefObj$numberOfTaxa() <
            ds@genotype@jRefObj$numberOfTaxa()
    )
    expect_equal(
        nrow(as.data.frame(sub)),
        sub@genotype@jRefObj$numberOfTaxa()
    )
})


# /// removeMinorSNPStates() /////////////////////////////////////////

test_that("removeMinorSNPStates() returns the input class", {
    expect_s4_class(removeMinorSNPStates(gt), "TasselGenotype")
    expect_s4_class(removeMinorSNPStates(ds), "TasselGenomicDataset")
})

test_that("removeMinorSNPStates() keeps taxa and sites intact", {
    out <- removeMinorSNPStates(gt)
    expect_equal(out@jRefObj$numberOfTaxa(), gt@jRefObj$numberOfTaxa())
    expect_true(out@jRefObj$numberOfSites() <= gt@jRefObj$numberOfSites())
})

test_that("removeMinorSNPStates() matches the legacy site filter", {
    withr::local_options(lifecycle_verbosity = "quiet")
    legacy <- filterGenotypeTableSites(
        readGenotypeTableFromPath(rtFiles$gt_hmp_path),
        removeMinorSNPStates = TRUE
    )
    legacyGt <- rTASSEL:::newTasselGenotype(getGenotypeTable(legacy), gt)
    expect_equal(as.matrix(removeMinorSNPStates(gt)), as.matrix(legacyGt))
})

test_that("removeMinorSNPStates() rejects unsupported input", {
    expect_error(removeMinorSNPStates(phenotype(ds)), "needs genotype data")
    expect_error(removeMinorSNPStates("not an rTASSEL object"), "Unsupported")
})


# /// Error handling /////////////////////////////////////////////////

test_that("bracket rejects invalid taxa selector types", {
    expect_error(gt[42L, ])
    expect_error(gt[TRUE, ])
})

test_that("bracket rejects invalid site selector types", {
    expect_error(gt[, TRUE])
})

test_that("region() errors for nonexistent chromosome", {
    expect_error(gt[, region("chrZZZ", 1, 100)])
})

test_that("taxa selection errors when no taxa match", {
    expect_error(gt[taxa("NONEXISTENT_TAXON_XYZ"), ])
})
