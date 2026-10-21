# === Tests for the dplyr-style filter verbs =========================
#
# The verbs are a second grammar over the machinery `[` already uses, so
# most of what is asserted here is that a verb and its bracket
# equivalent return the same taxa and sites.

## Preamble ----
startLogger()

# `lifecycle::expect_deprecated()` matches on condition class, which
# needs the third edition of testthat
local_edition(3)

gt <- rtObjs$gt_hmp
ds <- rtObjs$ds_hmp_ph_nomiss

nTaxa  <- function(x) getGenotypeTable(x)$numberOfTaxa()
nSites <- function(x) getGenotypeTable(x)$numberOfSites()

siteNamesOf <- function(x) rTASSEL:::batchSiteNames(getGenotypeTable(x))

testGr <- GenomicRanges::GRanges(
    seqnames = c("1", "2"),
    ranges   = IRanges::IRanges(
        start = c(1e6, 5e5),
        end   = c(5e6, 1e7)
    )
)


# /// filterSites() //////////////////////////////////////////////////

test_that("filterSites() matches the equivalent sitesWhere()", {
    verb <- gt |> filterSites(maf >= 0.05)
    expect_s4_class(verb, "TasselGenotype")
    expect_equal(siteNamesOf(verb), siteNamesOf(gt[, sitesWhere(maf >= 0.05)]))
})

test_that("filterSites() combines several predicates with &", {
    verb <- gt |> filterSites(maf >= 0.05, !isIndel)
    expect_equal(
        siteNamesOf(verb),
        siteNamesOf(gt[, sitesWhere(maf >= 0.05 & !isIndel)])
    )
})

test_that("filterSites() reproduces chrom() and region()", {
    expect_equal(
        siteNamesOf(gt |> filterSites(chrom %in% c("1", "5"))),
        siteNamesOf(gt[, chrom("1", "5")])
    )
    expect_equal(
        siteNamesOf(gt |> filterSites(chrom == "1", pos >= 1e6, pos <= 5e6)),
        siteNamesOf(gt[, region("1", 1e6, 5e6)])
    )
})

test_that("filterSites() agrees whether or not the plugin handles it", {
    # A literal threshold is pushed down to FilterSiteBuilderPlugin, the
    # same threshold held in a variable is not, so both paths are covered
    thresh <- 0.05
    shortCircuited <- gt |> filterSites(maf >= 0.05, het <= 0.1)
    fellBack       <- gt |> filterSites(maf >= thresh, het <= 0.1)

    expect_equal(
        siteNamesOf(shortCircuited),
        siteNamesOf(gt[, sitesWhere(maf >= 0.05 & het <= 0.1)])
    )
    expect_equal(siteNamesOf(fellBack), siteNamesOf(shortCircuited))
})

test_that("filterSites() combines predicates from different environments", {
    # Injecting a quosure gives the first predicate an environment the
    # second does not share, so each is evaluated on its own terms
    thresh <- 0.05
    viaWrapper <- function(x, predicate) {
        filterSites(x, !!predicate, het <= 0.1)
    }

    expect_equal(
        siteNamesOf(viaWrapper(gt, rlang::quo(maf >= thresh))),
        siteNamesOf(gt[, sitesWhere(maf >= 0.05 & het <= 0.1)])
    )
})

test_that("several predicates are still pushed down to TASSEL", {
    # Push-down is what keeps a verb call from pulling every site's
    # metadata into R, so assert the selector reaches the plugin rather
    # than inferring it from the result
    selector <- rTASSEL:::predicateSiteSelector(
        rlang::quos(maf >= 0.05, het <= 0.1)
    )
    pushed <- rTASSEL:::tryPluginShortCircuit(getGenotypeTable(gt), selector)

    expect_s4_class(pushed, "jobjRef")
    expect_equal(
        pushed$numberOfSites(),
        nSites(gt[, sitesWhere(maf >= 0.05 & het <= 0.1)])
    )
})

test_that("filterSites() rejects a non-logical predicate", {
    expect_error(gt |> filterSites(maf), "logical vector")
})

test_that("filterSites() errors when nothing is left", {
    expect_error(gt |> filterSites(maf >= 0.99), "No sites match")
})


# /// filterTaxa() ///////////////////////////////////////////////////

test_that("filterTaxa() matches the equivalent taxaWhere()", {
    verb <- gt |> filterTaxa(startsWith(taxaId, "CML"))
    expect_s4_class(verb, "TasselGenotype")
    expect_equal(
        taxaList(verb),
        taxaList(gt[taxaWhere(startsWith(taxaId, "CML")), ])
    )
})

test_that("filterTaxa() combines several predicates with &", {
    expect_equal(
        taxaList(gt |> filterTaxa(notMissing >= 0.8, het <= 0.1)),
        taxaList(gt[taxaWhere(notMissing >= 0.8 & het <= 0.1), ])
    )
})

test_that("filterTaxa() agrees whether or not the plugin handles it", {
    thresh <- 0.9
    shortCircuited <- gt |> filterTaxa(notMissing >= 0.9, het <= 0.01)
    fellBack       <- gt |> filterTaxa(notMissing >= thresh, het <= 0.01)

    expect_equal(taxaList(fellBack), taxaList(shortCircuited))
    expect_equal(
        taxaList(shortCircuited),
        taxaList(gt[taxaWhere(notMissing >= 0.9 & het <= 0.01), ])
    )
})

test_that("filterTaxa() errors when nothing is left", {
    expect_error(gt |> filterTaxa(notMissing >= 1.1), "No taxa match")
})


# /// selectSites() //////////////////////////////////////////////////

test_that("selectSites() matches siteIds() for literal names", {
    markers <- c("PZB00859.1", "PZA01271.1")
    expect_equal(
        siteNamesOf(gt |> selectSites("PZB00859.1", "PZA01271.1")),
        siteNamesOf(gt[, siteIds(markers)])
    )
})

test_that("selectSites() accepts tidyselect helpers", {
    markers <- c("PZB00859.1", "PZA01271.1")

    expect_equal(
        siteNamesOf(gt |> selectSites(all_of(markers))),
        siteNamesOf(gt[, siteIds(markers)])
    )
    expect_equal(
        siteNamesOf(gt |> selectSites(any_of(c(markers, "NOT_A_MARKER")))),
        siteNamesOf(gt[, siteIds(markers)])
    )

    prefixed <- gt |> selectSites(starts_with("PZA00"))
    expect_true(all(startsWith(siteNamesOf(prefixed), "PZA00")))
    expect_true(nSites(prefixed) > 0)
})

test_that("selectSites() drops with a leading minus", {
    markers <- c("PZB00859.1", "PZA01271.1")
    expect_equal(
        siteNamesOf(gt |> selectSites(-any_of(markers))),
        siteNamesOf(gt[, !siteIds(markers)])
    )
})

test_that("selectSites() errors on an absent all_of() ID", {
    expect_error(gt |> selectSites(all_of("NOT_A_MARKER")))
})

test_that("selectSites() errors when nothing is selected", {
    expect_error(gt |> selectSites(any_of("NOT_A_MARKER")), "No sites match")
})

test_that("selectSites() rejects renaming", {
    expect_error(gt |> selectSites(newName = "PZB00859.1"), "rename")
})


# /// selectTaxa() ///////////////////////////////////////////////////

test_that("selectTaxa() matches taxa() for literal IDs", {
    expect_equal(
        taxaList(gt |> selectTaxa("33-16", "38-11")),
        taxaList(gt[taxa("33-16", "38-11"), ])
    )
})

test_that("selectTaxa() accepts tidyselect helpers", {
    expect_equal(
        taxaList(gt |> selectTaxa(starts_with("CML"))),
        taxaList(gt[taxaWhere(startsWith(taxaId, "CML")), ])
    )
    expect_equal(
        taxaList(gt |> selectTaxa(matches("^CML"))),
        taxaList(gt[taxaWhere(grepl("^CML", taxaId)), ])
    )
})

test_that("selectTaxa() drops with a leading minus", {
    dropped <- c("33-16", "38-11")
    expect_equal(
        taxaList(gt |> selectTaxa(-any_of(dropped))),
        taxaList(gt[!taxa(dropped), ])
    )
})

test_that("selectTaxa() errors when nothing is selected", {
    expect_error(gt |> selectTaxa(any_of("NOT_A_TAXON")), "No taxa match")
})


# /// sliceSites() ///////////////////////////////////////////////////

test_that("sliceSites() matches sites() for positive positions", {
    expect_equal(
        siteNamesOf(gt |> sliceSites(1:100)),
        siteNamesOf(gt[, sites(1:100)])
    )
    expect_equal(
        siteNamesOf(gt |> sliceSites(c(1, 7, 42))),
        siteNamesOf(gt[, sites(c(1, 7, 42))])
    )
})

test_that("sliceSites() positions are 1-based", {
    expect_equal(siteNamesOf(gt |> sliceSites(1:10)), positionList(gt)$Name[1:10])
})

test_that("sliceSites() drops on negative positions", {
    expect_equal(
        siteNamesOf(gt |> sliceSites(-(1:10))),
        siteNamesOf(gt[, !sites(1:10)])
    )
})

test_that("sliceSites() takes positions across several arguments", {
    expect_equal(
        siteNamesOf(gt |> sliceSites(1:5, 10, 20)),
        siteNamesOf(gt[, sites(c(1:5, 10, 20))])
    )
})

test_that("sliceSites() ignores zeros and positions past the end", {
    total <- nSites(gt)

    expect_equal(
        siteNamesOf(gt |> sliceSites(c(0, 1, 2, total + 100))),
        siteNamesOf(gt[, sites(1:2)])
    )
    expect_equal(nSites(gt |> sliceSites(-(total + 1))), total)
})

test_that("sliceSites() cannot mix signs", {
    expect_error(gt |> sliceSites(c(1, -2)), "positive and negative")
})

test_that("sliceSites() rejects non-numeric positions", {
    expect_error(gt |> sliceSites("PZB00859.1"), "must be numeric")
    expect_error(gt |> sliceSites(1.5), "loss of precision")
})

test_that("sliceSites() errors when every position is past the end", {
    expect_error(gt |> sliceSites(nSites(gt) + 1), "selected nothing")
})


# /// sliceTaxa() ////////////////////////////////////////////////////

test_that("sliceTaxa() keeps taxa by position", {
    expect_equal(taxaList(gt |> sliceTaxa(1:10)), taxaList(gt)[1:10])
})

test_that("sliceTaxa() drops on negative positions", {
    expect_equal(
        taxaList(gt |> sliceTaxa(-(1:10))),
        taxaList(gt[!taxa(taxaList(gt)[1:10]), ])
    )
})

test_that("sliceTaxa() ignores positions past the end", {
    expect_equal(nTaxa(gt |> sliceTaxa(-(nTaxa(gt) + 1))), nTaxa(gt))
})


# /// overlaps() /////////////////////////////////////////////////////

test_that("overlaps() matches region() for the same ranges", {
    expect_equal(
        siteNamesOf(gt |> filterSites(overlaps(testGr))),
        siteNamesOf(gt[, region(testGr)])
    )
})

test_that("overlaps() works inside sitesWhere() too", {
    expect_equal(
        siteNamesOf(gt[, sitesWhere(overlaps(testGr))]),
        siteNamesOf(gt[, region(testGr)])
    )
})

test_that("overlaps() combines with other site criteria", {
    combined <- gt |> filterSites(overlaps(testGr), maf >= 0.05)
    expect_equal(
        siteNamesOf(combined),
        siteNamesOf(gt[, region(testGr)][, sitesWhere(maf >= 0.05)])
    )
})

test_that("overlaps() can be negated", {
    expect_equal(
        nSites(gt |> filterSites(!overlaps(testGr))),
        nSites(gt) - nSites(gt[, region(testGr)])
    )
})

test_that("overlaps() errors outside a site predicate", {
    expect_error(overlaps(testGr), "must be used inside")
})

test_that("overlaps() needs a GRanges object", {
    expect_error(gt |> filterSites(overlaps("chr1")), "<GRanges>")
})


# /// Chaining and no-ops ////////////////////////////////////////////

test_that("verbs chain to filter both axes", {
    piped <- gt |>
        filterTaxa(notMissing >= 0.8) |>
        filterSites(maf >= 0.05)
    bracketed <- gt[taxaWhere(notMissing >= 0.8), sitesWhere(maf >= 0.05)]

    expect_equal(taxaList(piped), taxaList(bracketed))
    expect_equal(siteNamesOf(piped), siteNamesOf(bracketed))
})

test_that("verbs called with no arguments return the input unchanged", {
    expect_identical(filterSites(gt), gt)
    expect_identical(filterTaxa(gt), gt)
    expect_identical(selectSites(gt), gt)
    expect_identical(selectTaxa(gt), gt)
    expect_identical(sliceSites(gt), gt)
    expect_identical(sliceTaxa(gt), gt)
})


# /// Return type preservation ///////////////////////////////////////

test_that("verbs return the class they were given", {
    numGt <- readGenotype(rtMatrices$num_gt_md)
    expect_s4_class(numGt |> sliceSites(1:5), "TasselNumericGenotype")
    expect_s4_class(numGt |> selectTaxa(starts_with("line_0")), "TasselNumericGenotype")
})

test_that("verbs re-join phenotype data on a genomic dataset", {
    sub <- ds |> filterTaxa(notMissing >= 0.9)

    expect_s4_class(sub, "TasselGenomicDataset")
    expect_true(nTaxa(sub) < nTaxa(ds))
    expect_equal(nrow(as.data.frame(sub)), nTaxa(sub))
})

test_that("verbs on a genomic dataset match the bracket form", {
    piped <- ds |>
        filterTaxa(notMissing >= 0.9) |>
        filterSites(maf >= 0.05)
    bracketed <- ds[taxaWhere(notMissing >= 0.9), sitesWhere(maf >= 0.05)]

    expect_equal(taxaList(piped), taxaList(bracketed))
    expect_equal(siteNamesOf(piped), siteNamesOf(bracketed))
})


# /// Input handling /////////////////////////////////////////////////

test_that("verbs reject input without genotype data", {
    expect_error(
        filterSites(rtObjs$ph_nomiss, maf >= 0.05),
        "needs genotype data"
    )
    expect_error(sliceTaxa("not an rTASSEL object", 1), "Unsupported")
})

test_that("verbs warn but work for TasselGenotypePhenotype input", {
    lifecycle::expect_deprecated(
        out <- filterSites(rtObjsLegacy$gt_hmp, maf >= 0.05),
        "TasselGenotypePhenotype"
    )
    expect_s4_class(out, "TasselGenotypePhenotype")
    expect_equal(
        siteNamesOf(out),
        siteNamesOf(gt[, sitesWhere(maf >= 0.05)])
    )
})
