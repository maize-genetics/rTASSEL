# === Tests for `seqDiversity()` ================================

## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasPheno   <- rtObjs$ph_nomiss
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss


## Error tests ----

### General errors
test_that("seqDiversity() throws general exceptions.", {
    expect_error(
        object = seqDiversity(mtcars),
        regexp = "Unsupported input object"
    )

    expect_error(
        object = seqDiversity(tasGeno, endSite = 5000000),
        regexp = "End site is out of bounds. Max index bound is: 3092"
    )
})


## Return tests ----
test_that("seqDiversity() returns correct data structures.", {
    seqDF <- seqDiversity(tasGeno)

    expect_equal(length(seqDF), 2)
    expect_equal(names(seqDF), c("Diversity", "PolyDist"))

    expect_equal(dim(seqDF$Diversity), c(1, 14))
    expect_equal(dim(seqDF$PolyDist), c(563, 2))

    expect_equal(
        object   = colnames(seqDF$Diversity),
        expected = c(
            "Site_Type", "Chromosome", "StartChrPosition", "EndChrPosition",
            "StartSite", "EndSite", "MidSite", "SiteCount", "AvgSiteCount",
            "SegSites", "PiPerBP", "ThetaPerBP", "Haplotypes", "TajimaD"
        )
    )

    expect_equal(
        object   = colnames(seqDF$PolyDist),
        expected = c("Site_Freq", "ALLs0-e3092")
    )
})


test_that("seqDiversity() accepts a genomic dataset", {
    seqDF <- seqDiversity(tasDataset)

    expect_equal(names(seqDF), c("Diversity", "PolyDist"))
    expect_equal(dim(seqDF$Diversity), c(1, 14))
})


## Back-compatibility ----
test_that("seqDiversity() accepts a deprecated TasselGenotypePhenotype", {
    expect_equal(
        seqDiversity(rtObjsLegacy$gt_hmp),
        seqDiversity(tasGeno)
    )
})



