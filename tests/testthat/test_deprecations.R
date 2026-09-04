# === Tests for the 0.14.0 deprecations ==============================
#
# `helper_vars.R` silences lifecycle notices for the rest of the suite,
# so every test here re-enables them - `lifecycle::expect_deprecated()`
# does that on its own, the remaining tests ask for it explicitly.

## Preamble ----
startLogger()

# `expect_deprecated()` matches on condition class, which needs the third
# edition of testthat
local_edition(3)

gt <- rtObjs$gt_hmp
ds <- rtObjs$ds_hmp_ph_nomiss
ph <- rtObjs$ph_nomiss

legacyGt  <- rtObjsLegacy$gt_hmp
legacyPh  <- rtObjsLegacy$ph_nomiss
legacyGtPh <- rtObjsLegacy$gt_hmp_ph_nomiss

# Two phenotypes carrying distinct traits, so the joins below have
# something to join on
smallPhDfs <- list(
    a = data.frame(taxa = c("a", "b", "c", "d"), weight = c(120, 150, 100, 70)),
    b = data.frame(taxa = c("a", "b", "c"), height = c(12, 15, 10))
)
smallPhAttr <- function(trait) {
    data.frame(
        col_id      = c("taxa", trait),
        tassel_attr = c("taxa", "data")
    )
}


# /// Deprecated readers /////////////////////////////////////////////

test_that("readGenotypeTableFromPath() points at readGenotype()", {
    lifecycle::expect_deprecated(
        readGenotypeTableFromPath(rtFiles$gt_hmp_path),
        "readGenotype\\(\\)"
    )
})

test_that("readPhenotypeFromPath() points at readPhenotype()", {
    lifecycle::expect_deprecated(
        readPhenotypeFromPath(rtFiles$ph_nomiss_path),
        "readPhenotype\\(\\)"
    )
})

test_that("readPhenotypeFromDataFrame() points at readPhenotype()", {
    lifecycle::expect_deprecated(
        readPhenotypeFromDataFrame(smallPhDfs$a, "taxa"),
        "readPhenotype\\(\\)"
    )
})

test_that("readGenotypePhenotype() points at readGenomicDataset()", {
    lifecycle::expect_deprecated(
        readGenotypePhenotype(rtFiles$gt_hmp_path, rtFiles$ph_nomiss_path),
        "readGenomicDataset\\(\\)"
    )
})

test_that("readGenotypePhenotype() reports itself, not its delegates", {
    withr::local_options(lifecycle_verbosity = "warning")

    warns <- capture_warnings(
        readGenotypePhenotype(rtFiles$gt_hmp_path, rtFiles$ph_nomiss_path)
    )

    expect_length(warns, 1L)
    expect_match(warns, "readGenotypePhenotype\\(\\)")
})


# /// Deprecated getters /////////////////////////////////////////////

test_that("getPhenotypeDF() points at as.data.frame()", {
    lifecycle::expect_deprecated(
        getPhenotypeDF(ph),
        "as.data.frame\\(\\)"
    )
})

test_that("getSumExpFromGenotypeTable() is deprecated", {
    lifecycle::expect_deprecated(
        getSumExpFromGenotypeTable(gt, verbose = FALSE),
        "getSumExpFromGenotypeTable\\(\\)"
    )
})


# /// Deprecated filters /////////////////////////////////////////////

test_that("filterGenotypeTableSites() points at bracket subsetting", {
    lifecycle::expect_deprecated(
        filterGenotypeTableSites(gt, siteMinAlleleFreq = 0.1),
        "sitesWhere"
    )
})

test_that("filterGenotypeTableTaxa() points at bracket subsetting", {
    lifecycle::expect_deprecated(
        filterGenotypeTableTaxa(gt, taxa = c("33-16", "38-11")),
        "taxaWhere"
    )
})

test_that("filterGenotypeTableBySiteName() points at bracket subsetting", {
    lifecycle::expect_deprecated(
        filterGenotypeTableBySiteName(gt, "PZB00859.1"),
        "siteIds"
    )
})

test_that("old filters warn once, about themselves, for legacy input", {
    withr::local_options(lifecycle_verbosity = "warning")

    warns <- capture_warnings(
        filterGenotypeTableSites(legacyGtPh, siteMinAlleleFreq = 0.1)
    )

    expect_length(warns, 1L)
    expect_match(warns, "filterGenotypeTableSites\\(\\)")
})


# /// Deprecated input class /////////////////////////////////////////

test_that("a TasselGenotypePhenotype input is deprecated, per function", {
    withr::local_options(lifecycle_verbosity = "warning")

    lifecycle::expect_deprecated(
        kinshipMatrix(legacyGtPh),
        "`kinshipMatrix\\(\\)`"
    )
    lifecycle::expect_deprecated(
        distanceMatrix(legacyGtPh),
        "`distanceMatrix\\(\\)`"
    )
    lifecycle::expect_deprecated(
        pca(legacyGt),
        "`pca\\(\\)`"
    )
})

test_that("the deprecation notice names the modern replacements", {
    lifecycle::expect_deprecated(
        kinshipMatrix(legacyGtPh),
        "readGenomicDataset"
    )
})

test_that("modern classes are accepted without any notice", {
    withr::local_options(lifecycle_verbosity = "warning")

    expect_length(capture_warnings(kinshipMatrix(gt)), 0L)
    expect_length(capture_warnings(kinshipMatrix(ds)), 0L)
    expect_length(capture_warnings(pca(gt)), 0L)
})

test_that("unsupported input reports the accepted classes", {
    expect_error(kinshipMatrix(mtcars), "Unsupported input object")
    expect_error(kinshipMatrix(mtcars), "<TasselGenomicDataset>")
    expect_error(kinshipMatrix(mtcars), "<TasselGenotypePhenotype>")
})

test_that("input missing the required component is reported, not warned", {
    expect_error(kinshipMatrix(ph), "needs genotype data")
    expect_error(assocModelFitter(gt, . ~ .), "needs phenotype data")
})


# /// Class in, class out ////////////////////////////////////////////

test_that("old filters return the class they were given", {
    withr::local_options(lifecycle_verbosity = "quiet")

    expect_s4_class(
        filterGenotypeTableSites(gt, siteMinAlleleFreq = 0.1),
        "TasselGenotype"
    )
    expect_s4_class(
        filterGenotypeTableSites(ds, siteMinAlleleFreq = 0.1),
        "TasselGenomicDataset"
    )
    expect_s4_class(
        filterGenotypeTableSites(legacyGtPh, siteMinAlleleFreq = 0.1),
        "TasselGenotypePhenotype"
    )

    expect_s4_class(
        filterGenotypeTableTaxa(gt, taxa = c("33-16", "38-11")),
        "TasselGenotype"
    )
    expect_s4_class(
        filterGenotypeTableTaxa(ds, taxa = c("33-16", "38-11")),
        "TasselGenomicDataset"
    )
    expect_s4_class(
        filterGenotypeTableTaxa(legacyGtPh, taxa = c("33-16", "38-11")),
        "TasselGenotypePhenotype"
    )

    expect_s4_class(
        filterGenotypeTableBySiteName(gt, "PZB00859.1"),
        "TasselGenotype"
    )
    expect_s4_class(
        filterGenotypeTableBySiteName(ds, "PZB00859.1"),
        "TasselGenomicDataset"
    )
    expect_s4_class(
        filterGenotypeTableBySiteName(legacyGtPh, "PZB00859.1"),
        "TasselGenotypePhenotype"
    )
})

test_that("old filters keep phenotype data attached to their input", {
    withr::local_options(lifecycle_verbosity = "quiet")

    sub <- filterGenotypeTableBySiteName(ds, c("PZB00859.1", "PZA01271.1"))

    expect_equal(sub@genotype@jRefObj$numberOfSites(), 2L)
    expect_equal(traitNames(sub), traitNames(ds))
})

test_that("imputation returns the class it was given", {
    withr::local_options(lifecycle_verbosity = "quiet")

    expect_s4_class(imputeNumeric(gt[, sites(1:10)]), "TasselGenotype")
    expect_s4_class(imputeNumeric(ds[, sites(1:10)]), "TasselGenomicDataset")
    expect_s4_class(
        imputeNumeric(
            filterGenotypeTableSites(
                legacyGtPh,
                siteRangeFilterType = "sites",
                startSite = 0,
                endSite = 10
            )
        ),
        "TasselGenotypePhenotype"
    )
})

test_that("mergeGenotypeTables() returns TGP only when every input was one", {
    withr::local_options(lifecycle_verbosity = "quiet")

    expect_s4_class(mergeGenotypeTables(list(gt, gt)), "TasselGenotype")
    expect_s4_class(
        mergeGenotypeTables(list(legacyGt, legacyGt)),
        "TasselGenotypePhenotype"
    )
    expect_s4_class(mergeGenotypeTables(list(gt, legacyGt)), "TasselGenotype")
})

test_that("joins return TGP only when every input was one", {
    withr::local_options(lifecycle_verbosity = "quiet")

    modernA <- readPhenotype(smallPhDfs$a, attr = smallPhAttr("weight"))
    modernB <- readPhenotype(smallPhDfs$b, attr = smallPhAttr("height"))
    legacyA <- readPhenotypeFromDataFrame(smallPhDfs$a, "taxa")
    legacyB <- readPhenotypeFromDataFrame(smallPhDfs$b, "taxa")

    expect_s4_class(intersectJoin(c(modernA, modernB)), "TasselPhenotype")
    expect_s4_class(
        intersectJoin(c(legacyA, legacyB)),
        "TasselGenotypePhenotype"
    )
    expect_s4_class(intersectJoin(c(modernA, legacyB)), "TasselPhenotype")

    expect_s4_class(unionJoin(c(modernA, modernB)), "TasselPhenotype")
    expect_s4_class(unionJoin(c(legacyA, legacyB)), "TasselGenotypePhenotype")
})


# /// Legacy pipelines still work ////////////////////////////////////

test_that("a fully legacy pipeline runs end to end", {
    withr::local_options(lifecycle_verbosity = "quiet")

    filtered <- filterGenotypeTableTaxa(
        legacyGtPh,
        taxa = head(getTaxaIDs(legacyGtPh), 30)
    )
    kin <- kinshipMatrix(filtered)
    res <- assocModelFitter(filtered, . ~ ., fitMarkers = TRUE, kinship = kin)

    expect_s4_class(filtered, "TasselGenotypePhenotype")
    expect_s4_class(kin, "TasselDistanceMatrix")
    expect_s4_class(res, "AssociationResults")
})

test_that("legacy and modern inputs give the same result", {
    withr::local_options(lifecycle_verbosity = "quiet")

    expect_equal(
        as.matrix(kinshipMatrix(ds)),
        as.matrix(kinshipMatrix(legacyGtPh))
    )
})

test_that("the deprecated class is still recognised by the getters", {
    expect_false(rJava::is.jnull(getGenotypeTable(legacyGtPh)))
    expect_false(rJava::is.jnull(getPhenotypeTable(legacyGtPh)))
    expect_false(rJava::is.jnull(getTaxaList(legacyGtPh)))
    expect_false(rJava::is.jnull(getPositionList(legacyGtPh)))
    expect_true(rJava::is.jnull(getPhenotypeTable(legacyGt)))
})
