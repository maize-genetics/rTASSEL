# === Tests for bracket-based filtering ==============================

## Preamble ----
startLogger()

gt <- readGenotype(rtFiles$gt_hmp_path)


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
    sel <- sites(0:99)
    expect_s4_class(sel, "SiteSelector")
    expect_equal(sel@type, "indices")
    expect_equal(sel@indices, 0:99L)
    expect_false(sel@negate)
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
    sel <- sites(0:9)
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


# /// Bracket filtering: sites ///////////////////////////////////////

test_that("[sites()] selects by 0-based index", {
    sub <- gt[, sites(0:99)]
    expect_s4_class(sub, "TasselGenotype")
    expect_equal(sub@jRefObj$numberOfSites(), 100L)
})

test_that("[numeric] selects sites by plain numeric vector", {
    sub <- gt[, 0:49]
    expect_s4_class(sub, "TasselGenotype")
    expect_equal(sub@jRefObj$numberOfSites(), 50L)
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
    sub <- gt[, !sites(0:9)]
    expect_equal(sub@jRefObj$numberOfSites(), nTotal - 10L)
})


# /// Combined taxa + sites //////////////////////////////////////////

test_that("bracket filters taxa and sites simultaneously", {
    sub <- gt[taxa("33-16", "38-11"), sites(0:49)]
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
    sub <- numGt[, sites(0:4)]
    expect_s4_class(sub, "TasselNumericGenotype")
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
