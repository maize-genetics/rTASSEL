# === Tests for `seqDiversity()` ================================

## Preamble - load data ----

### Start logging info
startLogger()

### Load hapmap data
genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)

### Read data - need only non missing data!
phenoPathFast <- system.file(
    "extdata",
    "mdp_traits_nomissing.txt",
    package = "rTASSEL"
)

### Create rTASSEL phenotype only object
tasPheno <- readPhenotypeFromPath(
    path = phenoPathFast
)

### Create rTASSEL genotype only object
tasGeno <- readGenotypeTableFromPath(
    path = genoPathHMP
)

### Create rTASSEL object - use prior TASSEL genotype object
tasGenoPhenoFast <- readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoPathFast
)


## Error tests ----

### General errors
test_that("seqDiversity() throws general exceptions.", {
    expect_error(
        object = seqDiversity(mtcars),
        regexp = "`tasObj` must be of class `TasselGenotypePhenotype`"
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



